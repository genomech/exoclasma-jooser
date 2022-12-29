from Bio.Seq import Seq
import bisect
import logging
import re
import subprocess

'''
Load fragment map as dictionary
'''
def LoadFragmentMap(RestrictionSitesMap):
	FragmentMap = {}
	with open(RestrictionSitesMap, 'rt') as MapFile:
		for Contig in MapFile:
			List = Contig[:-1].split(' ')
			FragmentMap[List[0]] = [ int(item) for item in List[1:] ]
	return FragmentMap

'''
A stream which sorts merged_nodups lines, compresses them and saves into file.
'''
def GetMergedNoDupsSorter(FileName):
	SortCommand = f'sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n | gzip -c > "{FileName}"'
	Stream = subprocess.Popen(SortCommand, shell = True, executable = 'bash', stdin = subprocess.PIPE)
	return Stream

'''
1) Bisect finds fragment position on its contig
2) All fragments count of previous contigs are summed
'''
def FindFragment(Chrom, Pos, FragmentMap):
	FragmentNumber = bisect.bisect(FragmentMap[Chrom], Pos) + 1
	assert FragmentNumber < len(FragmentMap[Chrom]), f'Read position goes out the contig'
	for Contig in FragmentMap.keys():
		if Contig == Chrom: break
		FragmentNumber += len(FragmentMap[Contig])
	return FragmentNumber

'''
Standard merged_nodups scheme:
strand1 chr1 pos1 fragment1 strand2 chr2 pos2 fragment2 mapq1 cigar1 seq1 mapq2 cigar2 seq2 name1 name2
'''
def CreateMergedNoDupsLine(Read1, Read2, Pos1, Pos2, FragmentMap):
	Strand = lambda is_reverse: 16 if is_reverse else 0
	Fragment = lambda Read, Pos
	Line = list()
	Line.append(str(Strand(Read1.is_reverse)))
	Line.append(str(Read1.reference_name))
	Line.append(str(Pos1))
	Line.append(('0' if FragmentMap is None else str(FindFragment(Read1.reference_name, Pos1, FragmentMap))))
	Line.append(str(Strand(Read2.is_reverse)))
	Line.append(str(Read2.reference_name))
	Line.append(str(Pos2))
	Line.append(('1' if FragmentMap is None else str(FindFragment(Read2.reference_name, Pos2, FragmentMap))))
	Line.append(str(Read1.mapping_quality))
	Line.append(str(Read1.cigarstring))
	Line.append(Read1.seq.__str__())
	Line.append(str(Read2.mapping_quality))
	Line.append(str(Read2.cigarstring))
	Line.append(Read2.seq.__str__())
	Line.append(str(Read1.query_name))
	Line.append(str(Read2.query_name))
	return (' '.join(Line) + '\n')

'''
1) Looking for every match of ligation site regexp in forward, then in reverse complement.
2) Taking their end points (relative to read start)
3) If end points are the same, it is a ligation point

I'm not sure every restrictase behaves like that (keeps its site untouched),
but each one, represented in exoclasma-index config file, does. At least for now.
If I'm wrong, please report here: https://github.com/regnveig/exoclasma-jooser/issues
'''
def FindLigation(Sequence, SiteRegExp):
	Forward = Sequence.__str__()
	Reverse = Sequence.reverse_complement().__str__()
	Length = len(Forward)
	ForwardEnds = [Match.span()[1] for Match in SiteRegExp.finditer(Forward)]
	ReverseEnds = [(Length - Match.span()[1]) for Match in SiteRegExp.finditer(Reverse)]
	LigationPoints = tuple([Item for Item in ForwardEnds if Item in ReverseEnds])
	return LigationPoints

'''
Extract:
- Read number
- Read primary or not
- Read contig
- Read position: its end if reverse-stranded, or its start otherwise.
The start and the end are adjusted taking soft and hard clips into consideration.
'''
def ReadTyping(Read, ChromSizes):
	Result = {}
	Result['Number'] = 1 if Read.is_read1 else 2
	Result['Primary'] = not (Read.is_secondary or Read.is_supplementary)
	Result['Chr'] = str(Read.reference_name)
	# Extract clipped ends
	SoftOrHardClipped = (4, 5)
	CigarFirst = item.cigar[0]
	CigarLast = item.cigar[-1]
	if not Read.is_reverse:
		# Handle start
		Start = Read.reference_start + 1
		if CigarFirst[0] in SoftOrHardClipped:
			Start -= CigarFirst[1]
			if Start <= 0: Start = 1
		Result['Pos'] = int(Start)
	else:
		# Handle end
		End = Read.reference_end
		if CigarLast[0] in SoftOrHardClipped:
			End += CigarLast[1]
			if End >= ChromSizes[item.reference_name]: End = ChromSizes[item.reference_name]
		Result['Pos'] = int(End)
	Result['Read'] = Read
	return Result

def CalculateDistance(Coord1, Coord2):
	return float('+inf') if (Item1['Chr'] != Item2['Chr']) else abs(Item1['Pos'] - Item2['Pos'])

def SortItems(Item1, Item2):
	return tuple([(item['ID'], item['Pos']) for item in sorted([Item1, Item2], key = lambda x: (x['RefID'], x['Pos']))])

def GetDuplicationTag(Record):
	if not Record.is_duplicate: return False
	return dict(Record.tags)['DT']

def ProcessQuery(Query, ChromSizes):
	# Filter & count duplicates
	DuplicationTags = [GetDuplicationTag(Item) for Item in Query['ReadBlock']]
	if any([Item == 'SQ' for Item in DuplicationTags]):
		Query['Type'] = 'Optical Duplicates'
		return Query
	if any([Item == 'LB' for Item in DuplicationTags]):
		Query['Type'] = 'PCR Duplicates'
		return Query
	# Filter Unmapped
	if any([Item.is_unmapped for Item in Query['ReadBlock']]):
		Query['Type'] = 'Unmapped'
		return Query
	# Create Sorter
	TypeDict = { index: list() for index in ((1, True), '1s', '2p', '2s') }
	# Annotation
	for index, item in Query['ReadBlock']:
		Start = item.reference_start + 1
		End = item.reference_end
		CigarFirst = item.cigar[0]
		CigarLast = item.cigar[-1]
		SoftHard = (4, 5)
		if CigarFirst[0] in SoftHard:
			Start -= CigarFirst[1]
			if Start <= 0: Start = 1
		if CigarLast[0] in SoftHard:
			End += CigarLast[1]
			if End >= ChromSizes[item.reference_name]: End = ChromSizes[item.reference_name]
		Type = ('1' if item.is_read1 else '2') + ('s' if (item.is_secondary or item.is_supplementary) else 'p')
		TypeDict[Type].append({ 'ID': int(index), 'Chr': str(item.reference_name), 'RefID': int(item.reference_id), 'Pos': int(End) if item.is_reverse else int(Start) })
	# Create Pattern
	Pattern = tuple([len(item) for index, item in TypeDict.items()])
	TypeDict = { index: (None if not item else (item[0] if len(item) == 1 else item)) for index, item in TypeDict.items() }
	Dist = { f'1{index1}2{index2}': CalculateDistance(TypeDict[f'1{index1}'], TypeDict[f'2{index2}']) for index1, index2 in ('pp', 'ps', 'sp', 'ss')}
	# Norm Chimera 4 Ends
	if Pattern == (1, 1, 1, 1):
		if ((Dist['1p2p'] < 1000) and (Dist['1s2s'] < 1000)) or ((Dist['1p2s'] < 1000) and (Dist['1s2p'] < 1000)):
			Sorted = SortItems(TypeDict['1p'], TypeDict['1s'])
			Pair = [{ 'Read': Query['ReadBlock'][Sorted[0][0]][1], 'Pos': Sorted[0][1] }, { 'Read': Query['ReadBlock'][Sorted[1][0]][1], 'Pos': Sorted[1][1] }]
			return { 'ReadBlock': Query['ReadBlock'], 'Type': 'ChimericPaired', 'Pair': Pair }
		else: return { 'ReadBlock': Query['ReadBlock'], 'Type': 'ChimericAmbiguous' }
	# Norm Chimera 3 Ends
	elif Pattern in ((1, 0, 1, 1), (1, 1, 1, 0)):
		if TypeDict['1s'] is None:
			if ((Dist['1p2p'] < 1000) or (Dist['1p2s'] < 1000)): Sorted = SortItems(TypeDict['1p'], TypeDict['2p'] if Dist['1p2p'] > Dist['1p2s'] else TypeDict['2s'])
			else: Sorted = None
		if TypeDict['2s'] is None:
			if ((Dist['1p2p'] < 1000) or (Dist['1s2p'] < 1000)): Sorted = SortItems(TypeDict['2p'], TypeDict['1p'] if Dist['1p2p'] > Dist['1s2p'] else TypeDict['1s'])
			else: Sorted = None
		if Sorted is None: return { 'ReadBlock': Query['ReadBlock'], 'Type': 'ChimericAmbiguous' }
		Pair = [{ 'Read': Query['ReadBlock'][Sorted[0][0]][1], 'Pos': Sorted[0][1] }, { 'Read': Query['ReadBlock'][Sorted[1][0]][1], 'Pos': Sorted[1][1] }]
		return { 'ReadBlock': Query['ReadBlock'], 'Type': 'ChimericPaired', 'Pair': Pair }
	# Regular Pair
	elif Pattern == (1, 0, 1, 0):
		Sorted = SortItems(TypeDict['1p'], TypeDict['2p'])
		Pair = [{ 'Read': Query['ReadBlock'][Sorted[0][0]][1], 'Pos': Sorted[0][1] }, { 'Read': Query['ReadBlock'][Sorted[1][0]][1], 'Pos': Sorted[1][1] }]
		return { 'ReadBlock': Query['ReadBlock'], 'Type': 'NormalPaired', 'Pair': Pair }
	# Collisions
	elif (Pattern[1] > 1) or (Pattern[3] > 1):
		pass # TODO Collisions
	# Other
	return { 'ReadBlock': Query['ReadBlock'], 'Type': 'ChimericAmbiguous' }

def JooserFunc(**kwargs):
	N = types.SimpleNamespace(**kwargs)
	Input = pysam.AlignmentFile(N.Input_BAM, 'r', check_sq=False)
	Output = GetMergeNoDupsSorter(N.MergedNoDups_File)
	if N.Restriction_Site_Map is not None: FragmentMap = LoadFragmentMap(N.Restriction_Site_Map)
	ChromSizes = { Input.references[i]: Input.lengths[i] for i in range(Input.nreferences) }
	if N.Restriction_Site is not None: LigationQuery = re.compile(CONFIG_RESTRICTION_ENZYMES['Ligation'][N.Restriction_Site])
	Stats = { 'SequencedReadPairs': 0, 'NormalPaired': 0, 'ChimericPaired': 0, 'ChimericAmbiguous': 0, 'MappingQualityFailed': 0, 'PcrDuplicates': 0, 'OpticalDuplicates': 0, 'Unmapped': 0, 'Ligation': { 'Motif': None, 'LineCount': 0, 'PresentCount': 0 } }
	Query = { 'ReadName': None, 'ReadBlock': [] }
	def BlockProcess():
		Stats['SequencedReadPairs'] += 1
		#Stats["Ligation"]["LineCount"] += 1
				#if N.Restriction_Site is not None:
					#if LigationQuery.match() # TODO
		Query['ReadBlock'] = list(enumerate(Query['ReadBlock']))
		Result = ProcessQuery(Query, ChromSizes, N.Min_MAPQ)
		Stats[Result['Type']] += 1
		if Result['Type'] in ('ChimericPaired', 'NormalPaired'):
			Read1, Read2 = Result['Pair']
			Line = ' '.join([
				'16' if Read1['Read'].is_reverse else '0',
				str(Read1['Read'].reference_name),
				str(Read1['Pos']),
				'0' if N.Restriction_Site_Map is None else str(bisect.bisect(FragmentMap[Read1['Read'].reference_name], Read1['Pos'])),
				'16' if Read2['Read'].is_reverse else '0',
				str(Read2['Read'].reference_name),
				str(Read2['Pos']),
				'1' if N.Restriction_Site_Map is None else str(bisect.bisect(FragmentMap[Read2['Read'].reference_name], Read2['Pos'])),
				str(Read1['Read'].mapping_quality),
				str(Read1['Read'].cigarstring),
				str(Read1['Read'].seq.__str__()),
				str(Read2['Read'].mapping_quality),
				str(Read2['Read'].cigarstring),
				str(Read2['Read'].seq.__str__()),
				str(Read1['Read'].query_name),
				str(Read2['Read'].query_name)
				]) + '\n'
			Output.stdin.write(Line.encode('utf-8'))
	logging.info(RenderParameters('Start processing', N.Input_BAM))
	while 1:
		try:
			Record = next(Input)
			if Record.query_name == Query['ReadName']: Query['ReadBlock'].append(Record)
			else:
				BlockProcess()
				Query['ReadName'] = Record.query_name
				Query['ReadBlock'].clear()
				Query['ReadBlock'].append(Record)
		except StopIteration:
			BlockProcess()
			Input.close()
			Output.stdin.close()
			Output.wait()
			Stats['Alignable'] = Stats['ChimericPaired'] + Stats['NormalPaired']
			Stats['Duplicates'] = Stats['PcrDuplicates'] + Stats['OpticalDuplicates']
			for stat in ('ChimericPaired', 'ChimericAmbiguous', 'NormalPaired', 'Unmapped', 'Alignable', 'MappingQualityFailed', 'Duplicates', 'PcrDuplicates', 'OpticalDuplicates'): Stats[stat] = { 'Count': Stats[stat], '%': Stats[stat] / Stats['SequencedReadPairs'] * 100 }
			Stats['Ligation']['%'] = Stats['Ligation']['PresentCount'] / Stats['SequencedReadPairs'] * 100 # BUG WTF?
			# TODO Postprocessing? Library Complexity?
			json.dump(Stats, open(N.Jooser_Stats, 'wt'), indent=4, ensure_ascii=False)
			logging.info(RenderParameters('End processing', N.Input_BAM))
			break

def JuicerTools(**kwargs):
	logging.info(RenderParameters('*', '*'))
	for Key, Value in kwargs.items(): logging.info(RenderParameters(Key, Value))
	N = types.SimpleNamespace(**kwargs)
	RSString = '' if N.Restriction_Site_Map is None else f'-f "{N.Restriction_Site_Map}"'
	Command = f'java -jar "{JUICERTOOLS_PATH}" pre -j {N.Threads} {RSString} "{N.MergedNoDups_File}" "{N.Output_HIC_File}" "{N.ChromSizes_File}"'
	BashSubprocess(Name = 'JuicerTools', Command = Command)

def main():
	pass

if __name__ == '__main__': main()
