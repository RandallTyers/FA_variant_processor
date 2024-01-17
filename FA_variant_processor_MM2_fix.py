#!/usr/bin/env python3

"""
Start with using vast to analyze colony NGS data. 
	Use the F-sequence as the reference and the A-sequence as the cassette.
	Run this script using python3.

Goal: Analyze data from vast output and preserve strand connectivity while counting sequences. 
		Uses output from `vast report variant_frequency --groupby +`.

Steps:	
		1. Load ouptut from `vast report variant_frequency --groupby [tag(s) | +{=all}]
		2. 




# Usage: ./FA_variant_processsor.py 
# 
# 	e.g. python FA_variant_processsor.py "Elena_db/Reports/2021-12-22_21-03-50_variant_frequency_groupby_Library.Genotype.Tailed.Colony/F_Elena-REF-P.frequency.tsv"
# 
# 
# Outputs: 																							# TODO: Update — this is for FA_processor.py!
# 		1. Merged TSV files of Templated SNP frequencies (only data containing columns)
# 			columns: positions -- rows: file names
# 		2. TSV file of inferred genotypes based on CUFFOFF values (below)
# 			columns: positions -- rows: file names
# 		3. TSV file of inferred genotype frequencies
# 			columns: positions -- rows: inferred genotypes	
"""

'''
TODOs: 

I) Display information by position (for comparison with Elena's usual display)	√

II) Output tables in an 'attractive' format.

III) Bring histograms from FA_processor.py into this version.

# I) Display connectivity 
# 	• ¿simply pairs of colored lines for reference and alternate (cassette) alleles
# 
# II) Statistics
# 	• ¿normality e.g. when calling F|A?
# 	• ¿model in a way that accounts for contamination?
# 	• ¿use non-templated SNP data to estimate errors — expect 'similar' errors at templated sites?
# 
# III) Make reports optional? 
	• 
'''
####################################################################################################
##### Preset values ################################################################################
####################################################################################################

# 
INDEX_NAMES = ['Variant', 'Number of switches'] # TODO: why is read_csv (&ct) not preserving this?

# built from INDEX_NAMES —— last line in header of tsv output by vast report variant_frequency
LAST_HEADER_LINE = '\t'.join(INDEX_NAMES)	#"Variant\tNumber of switches"

EXPECTED_SITE_COUNT = 7

GENOTYPES = [
				"F",	# reference allele
				"F|A?",	# might have been heteroduplex or reference allele
				"F|A",	# heteroduplex
				"F?|A",	# might have been heteroduplex or cassette allele
				"A"		# cassette allele
]

# Cutoff values for discriminating among F F/A A
#	these are totally arbitrary, would be ideal to have some known data to get an idea of what works
F_CUTOFF = 0.15
FA_LOW_CUTOFF = 0.35
FA_HIGH_CUTOFF = 0.65
A_CUTOFF = 0.85

# alternative
TOLERANCE = 0.15

# F_CUTOFF = TOLERANCE
# FA_LOW_CUTOFF = 0.5 - TOLERANCE
# FA_HIGH_CUTOFF = 0.5 + TOLERANCE
# A_CUTOFF = 1.0 - TOLERANCE

CORR_CUTOFF = 0.5	# no longer in use?
KEEP_CUTOFF = 0.19	#0.13 #0.216	# 0.14 to 0.24 ...

CUTOFF_LABEL = f"KEEP_CUTOFF={KEEP_CUTOFF}"														# TODO: set this to "" to eliminate this once a cutoff value chosen?

## get top sequences —— currently all those with at least 1 occurence of a proportion > 0.05 ##
# minimum proportion to count
PROPORT_CUTOFF = 0.05	# was 0.1
# minimum counts to include
COUNT_CUTOFF = 1

# string used to indicate separation between strands in heteroduplex sequences
STRAND_SEP = ' ' * 3

# names (for levels of column headings) used to categorize data —— 'Library' dispensable later? 	# TODO: automate this ...
LEVELS = [
#			'Version',
#			'Library',
			'Genotype',
			'Tailed'
]

QUARTILES = (1,0.75,0.5,0.25,0)

### formating presets ###
SPACER_HEAVY = "=" * 90
SPACER_LIGHT = "-" * 90


####################################################################################################
##### Imports ######################################################################################
####################################################################################################

import numpy
import pandas
import matplotlib.pyplot as plt

import sys
import re
from os.path import isfile
from operator import itemgetter



# globally set plotting paramaters (TODO: figure out how to just change the x-axis tick labels)
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.size'] = '12'


	
####################################################################################################
##### Function definitions #########################################################################
####################################################################################################

# fix format of vast variant_frequency report tsv so it is readable
#	problem is that there are weird gaps (e.g. main header row only contains two labels with nothing
#		for all the data columns! (WTH were they thinking ...) 
def header_line_number(path, head_start=LAST_HEADER_LINE):
	"""
	Scans through lines of 'path' and returns the line number of first line matching 'head_start'.
	
	Source: https://stackoverflow.com/a/53636453/3613297 by d_kennetz.
	"""
	with open(path, "r") as fh:
		for line_num, line in enumerate(fh):
			if line.startswith(head_start):
				return line_num



# first fix format of vast variant_frequency report tsv so it is readable by pandas.read_csv
#	problem is that there are weird gaps (e.g. main header row only contains two labels with nothing
#		for all the data columns! (WTH were they thinking ...) 
def fixed_variants_and_header_len(variants_file):
	"""
	Creates a 'fixed' version of vast report variant_frequency if not already present. 
		Needed to make vast output parsable by pandas.read_csv.
	
	Input: Name of output from vast report variant_frequency.
	
	Output: Tuple contianing name of 'fixed' file and the length of the header.
	
	Requires: header_line_number, isfile (from os.path), re
	"""
	# if argument is a file name that looks like it was made by this function 
	#	then use argument as the "fixed" name (and don't try to fix it) 
	if "_fixed." in variants_file:
		fixed_variants_file = variants_file
		header_len = header_line_number(fixed_variants_file) + 1
	
	# if file name (argument) doesn't look "fixed" then create "fixed" file name	
	else:
		fullname, suffix = variants_file.rsplit(sep=".", maxsplit=1)
		
		fixed_variants_file = fullname + "_fixed." + suffix
		
		# check that this hasn't already been done (could happen if unfixed name used by mistake)
		if isfile(fixed_variants_file):
			# if already "fixed", just figure out the header length
			header_len = header_line_number(fixed_variants_file) + 1						 		# now getting: TypeError: unsupported operand type(s) for +: 'NoneType' and 'int'
			
		else: 
			# make a "fixed" version
			with open(variants_file, "r") as in_fh, open(fixed_variants_file, "w") as out_fh:
				
				## REGEX pattern for extracting header from VAST variant_frequency report ##
				#	probably should be using raw strings here ...
				#
				# Match header label information #
				#	groups: tag, value(s)
				#	where:
				#		tag is a data category (see VAST for usage)
				#		value(s) are the possibilities for that tag (as TSVs)
				header_pattern = re.compile("([^\t]+)\t\t(.*)")
				
				# parse and separate header from body
				#	group(1) includes header (or blank) lines up to and including 'LAST_HEADER_LINE'
				#	group(2) includes everything after the header
				match = re.fullmatch(
										f"((?:[^\t]+\t\t[^\n]*\n|\n)+{LAST_HEADER_LINE}\s*\n)(.*)", 
										in_fh.read(),
										flags=re.DOTALL
									)
				
				## parse header information into a dictionary ##
				# dictionary stores information for header labels
				#	keys: tags
				#	values: lists of tag value(s)
				heads_dict = dict([
									( head_match.group(1).strip(), head_match.group(2).split() ) 
										for head_match in header_pattern.finditer(match.group(1) ) 	# (SOLVED ASCII vs. utf-8 & ∆ superfluous tabs) was getting: AttributeError: 'NoneType' object has no attribute 'group'
				])
				
				# write tags & values into "fixed" version of tsv file
				#	TODO: ¿more efficient to go directly to dataframe rather than writing to disc?
				for tag, values in heads_dict.items():
					
					tsv_string = "\t".join( str(val) for val in values )
					
					out_fh.write( f"{tag}\t\t{tsv_string}\n" )
				
				# number of lines in header
				header_len = len(heads_dict) + 1
				
				# remove 'Number of Reads' data so it isn't in "condensed" headers
				del heads_dict['Number of Reads']
				
				## make last (main) header line 
				#	each column gets a "condensed" header —— "tag1value tag2value ... tag'n'value"
				# 		• '*heads_dict.values()' unpacks dictionary into lists
				#		• zip then groups across the lists
				#		• each group is joined into a string separated by spaces
				#		• those strings are then joined into a string separated by tabs
				header_tail = "\t".join( 
											" ".join( str(item) for item in items ) 
											for items in zip(*heads_dict.values()) 
				)
				
				# write condensed tags into "fixed" version of tsv file == main header line
				out_fh.write( f"{LAST_HEADER_LINE}\t{header_tail}\n" )
				
				# dump remainder of file (the data) into the new "fixed" tsv file
				out_fh.write( match.group(2) )
				
	
	return (fixed_variants_file, header_len)



# remove non-sense results due to vast errors #
# type 1: multiple tracts — signatures: 
#										more than one tuple of tuples inside the frozen set 
#										multiple (comma separated) values for "Number of switches"
def remove_vast_errors(df):
	"""
	Eliminate and report problem records from vast output.

	Currently recognizing errors where mutliple tracts appear to be found, these include finding:
		more than one tuple of tuples inside the frozen set
		multiple (comma separated) values for "Number of switches"
	
	Input is a dataframe with two level index where:
		the first level is a frozenset, and
		the second level should be an integer. 
	
	
	"""
	detailed_errors = []
	
	def is_numeric(object):
		"""
		Returns whether passed object is or can be converted into an int.
		"""
		result = False
		
		try:
			_ = int(object)
		except (TypeError, ValueError) as e:
			detailed_errors.append(e)
		else:
			result = True
		
		return result
	
	
	vast_errors = [
					fset for fset, switch_num in df.index
							if len(eval(fset)) > 1 or not is_numeric(switch_num)					# problem was that most (but not all) of switch_num types are str
	]
	
	if vast_errors:
		vast_error_out = f"vast_errors_{out_file_base}"
		print(
				f"Removing {len(vast_errors)} rows from data that appear to be vast errors.\n"
				f"Please check {vast_error_out} to make sure that this was the correct decision.",
				file=sys.stderr
		)
		write_df_to_tsv(df.loc[vast_errors], vast_error_out, index_label=None)
		
		with open(f"detailed_{vast_error_out}", mode='w') as error_out:								# was mode='x'
			error_out.write('\n'.join(repr(err) for err in detailed_errors))
		
		df.drop(labels=vast_errors, axis='index', level=0, inplace=True)
	
	return df



def write_df_to_tsv(dataframe, out_file_name, index_label='File name'):
	"""
	Wrapper for pandas.DataFrame.to_csv() to create tsv formatted files.
	
	Warns (with opportunity to abort) about overwriting files.
	
	Prints message to stderr after saving completed.
	
	Requires: sys, isfile, pandas
	"""
	# warn/abort if out_file exists
	if isfile(out_file_name):
		try:
			assert( 
					input(f"WARNING: '{out_file_name}' already exists — overwrite? (Y/N)").lower() 
					in ("yes", "y") 
			)
		except AssertionError:
			sys.exit("Aborting execution — please rename output files you wish to keep.")
	
	# write dataframe contents to file as tsv
	dataframe.to_csv( 
						path_or_buf=out_file_name, 
						sep='\t', 
						index_label=index_label
	)
	print(f"saved `{out_file_name}`", file=sys.stderr)



def save_show_plot(out_file_name):
	"""
	Wrapper for matplotlib.pyplot.savefig().
	
	Warns (with opportunity to abort) about overwriting files.
	
	Prints message to stderr after saving completed.
	
	Finally shows figure.
	
	Requires: sys, isfile, pandas, matplotlib.pyplot as plt
	"""
	# warn/abort if out_file exists
	if isfile(out_file_name):
		try:
			assert( 
					input(f"WARNING: '{out_file_name}' already exists — overwrite? (Y/N)").lower() 
					in ("yes", "y") 
			)
		except AssertionError:
			sys.exit("Aborting execution — please rename output files you wish to keep.")
	
	plt.savefig(out_file_name)
	
	print(f"saved `{out_file_name}`\n", file=sys.stderr)
	
	# show message
# 	print("Showing saved figure — close the figure window to continue.")
# 	
# 	plt.show()
	plt.close()


 
# TODO: Generalize√, add error checking i.e. seqs are same length (and at most two?)
def decompose(seq_string, strand_sep=STRAND_SEP, sep_char=None):
	"""
	Separates a 'best2' sequence into components.
	
	Input: potentially heteroduplex sequences as a string.
		 	Individual characters are separated by 'sep_char' — default (None) uses any white space.
		 	Separate strands are separated by 'strand_sep' — default is three spaces ('   ').
	
	Output: tuple of frozensets where each frozenset contains one or two characters.
			Position of frozenset within tuple corresponds to position within sequence.
	"""
	return	tuple(
					frozenset(chars)
					for chars in zip( *[
											seq.split(sep=sep_char)
											for seq in seq_string.split(sep=STRAND_SEP)
										]
					)
	)


# def xs_gt_by_pos_at_Allele(key):																	# <<<<<<<< mysteriously this throws a KeyError the first time the script is run, but then works if rerun (happens everytime a new version is saved ...)
# 	"""
# 	Wrapper for xs applied to genotype_by_position_df.
# 	
# 	Input: key — value of level 'Allele'.
# 	
# 	Output: column "slice" of genotype_by_position_df where level=='Allele' & value is key.
# 	
# 	Requires: pandas.
# 	"""
# 	return genotype_by_position_df.xs(
# 										key=key,
# 										axis='columns',
# 										level='Allele',
# 	) 

# # TODO: Fix this, it isn't picking best two sequences ...
# def best_indices_to_header(	
# 							df, 
# 							row_corr_df=pandas.DataFrame(), 
# 							CORR_CUTOFF=CORR_CUTOFF, 
# 							KEEP_CUTOFF=KEEP_CUTOFF
# ):
# 	"""
# 	Inputs: dataframe and dataframe of the inter row correlations for that dataframe plus cutoffs.
# 			• row_corr_df: correlations between rows, automatically calculated if not given
# 			• CORR_CUTOFF: highest correlation that won't be used to merge rows
# 			• KEEP_CUTOFF: hightest value (genotype frequency) that won't be included in header
# 	
# 	Returns: a collapsed version of the dataframe that has has a level added to the header. 
# 	
# 	Requires: pandas, numpy, 
# 	"""
# 	#CORR_CUTOFF = 0.5
# 	#KEEP_CUTOFF=0.14
# 	df = df.copy()
# 	
# 	if row_corr_df.empty:
# 		row_corr_df = df.T.corr()
# 	
# 	# analyze each column (colony) separately ––– current version is quite slow ... TODO: speed up?
# 	for col in df.columns:
# 		
# 		# look at each row (sequence variant) starting with most common 
# 		for row in df[col].sort_values(ascending=False).index:
# 			
# 			# get top two largest (largest will be self) correlations as a pandas Series
# 			#two_largest_corrs = row_corr_df[row].nlargest(2, row)
# 			# above fails, but below works —— bug in Series implementation of nlargest	?????
# 			two_largest_corrs = row_corr_df.nlargest(2, row)[row]
# 		
# 			# is this correlation better than the (arbitrary) cutoff value?
# 			if two_largest_corrs[1] > CORR_CUTOFF:
# 				# most correlated other row (sequence variant), putative PCR induced variant of row
# 				corr_row = two_largest_corrs.index[1]
# 			
# 				# is the putative PCR variant at least as complex as its supposed parent?
# 				if corr_row in df.index:
# 					# check to make sure adding less frequent to more frequent sequence
# 					if df.loc[row, col] < df.loc[corr_row, col]:									#(line 39) — problems if there NaN values?
# 						# switch row and corr_row so that sequence with larger proportion is being kept
# 						row, corr_row = corr_row, row
# 					# check that putative PCR variant at least as complex as its supposed parent
# 					if corr_row[1] >= row[1]: 	
# 						# add value of variant to putative parent
# 						df.loc[row, col] += df.loc[corr_row, col]	#[0]							#(line 45) — problems if there NaN values?
# 						# replace variant proportion with NaN
# 						df.loc[corr_row, col] = numpy.NaN
# 	
# 	# purge all rows that now entirely NaN
# 	df.dropna(how='all', inplace=True)
# 	
# 	head_df = df.columns.to_frame()
# 	head_df['best2'] = [
# 						tuple(
# 								sorted(		# sort simplifies results, but loses relative frequency
# 										(
# 											tup[0] for tup in df.nlargest(2, col).index.to_list()
# 													if df[col].loc[tup] > KEEP_CUTOFF				# weirdness going on here: more rational versions df.loc[tup, col] & df[col][tup] don't work!
# 										),
# 										reverse=True
# 								)
# 						) 
# 						for col in df.columns
# 	]
# 	
# 	df.columns = pandas.MultiIndex.from_frame(head_df)
# 	
# 	return df


####################################################################################################
##### MAIN #########################################################################################
####################################################################################################

# doubly sorted to deal with numbered files that have unequal numbers of digits (i.e. lengths)
variants_file = sys.argv[1]			#variants_file = "Elena_db/Reports/2021-12-22_21-03-50_variant_frequency_groupby_Library.Genotype.Tailed.Colony/F_Elena-REF-P.frequency.tsv"

	#variants_file = "Elena_db/Reports/2022-01-09_20-00-37_variant_frequency_groupby_Library.Genotype.Tailed.Colony_exclude_Library_WTmsh2MSH6tail_S1/F_Elena-REF-P.frequency.tsv"
	#variants_file = "Elena_db/Reports/2022-01-20_18-44-04_variant_frequency_groupby_Library.Genotype.Tailed.Colony_exclude_Library_WTmsh2MSH6Tailless_S2.WTmsh2MSH6tail_S1/F_Elena-REF-P.frequency.tsv"
	#variants_file = "Elena_db/Reports/2022-03-10_00-00-46_variant_frequency_groupby_Library.Genotype.Tailed.Colony_include_Version_xsL280q35u10/F_Elena-REF-P.frequency.tsv"
	#variants_file = "Elena_db/Reports/2022-03-08_14-15-06_variant_frequency_groupby_Library.Genotype.Tailed.Colony/F_Elena-REF-P.frequency.tsv"
out_file_base = sys.argv[0].rpartition('/')[2] + '_on_' + variants_file.rpartition('/')[2]


### import desired data from tsv files #############################################################

# get header_len & fix the vast variant_frequency report (if not already present)
# first fix format of vast variant_frequency report tsv so it is readable by pandas.read_csv
#	problem is that there are weird gaps (e.g. main header row only contains two labels with nothing
#		for all the data columns! (WTH were they thinking ...) 
fixed_variants_file, header_len = fixed_variants_and_header_len(variants_file)

# make a dataframe containing only the desired row of data from each tsv file
#	TODO: figure out how to get the indexing done correctly in one step ... 
variants_df = pandas.read_csv( 
								fixed_variants_file, 
								sep='\t', 
								index_col=[0,1], 
								header=[*range(header_len)] 
							)

# add in missing names for levels of MultiIndex —— ¿haven't found a way to make read_csv do this ...
variants_df.index.names = INDEX_NAMES

### process data from tsv file #####################################################################

# limit dataframe to only those rows with data —— possibly superfluous (aka. purged variant dataframe)
purg_var_df = variants_df.loc[variants_df.max(1) > 0]	#.drop_duplicates() <- ¿mistake?

# sort rows by max row value
#purg_var_df = purg_var_df.reindex(purg_var_df.max(1).sort_values(ascending=False).index, axis='index')

# sort rows by number of columns containing a proportion greater than 0.1
purg_var_df = purg_var_df.reindex(
										purg_var_df[purg_var_df > 0.1]
											.count(axis='columns')
											.sort_values(ascending=False)
											.index,
										axis='index'
)

purg_var_df = remove_vast_errors(purg_var_df)


# # deal with duplicate indices (errors produced by vast) --- found nothing!!!! #
# duplicate_df = purg_var_df.iloc[purg_var_df.index.duplicated(keep=False)]
# 
# # check if there are any rows in this dataframe
# if len(duplicate_df.index) != 0:
# 	# if so, warn and save them to a tsv file for examination
# 	duplicate_out = "duplicates_" + out_file_base
# 	print(
# 			"Duplicate rows (sequences) found in the data, all, but the first (top) were deleted.\n"
# 			f"Please check {duplicate_out} to make sure that this was the correct decision."
# 	)
# 	write_df_to_tsv(duplicate_df, duplicate_out, index_label=None)
# 	
# 	# do the removal by keeping everything except the later duplicates 
# 	purg_var_df = purg_var_df[~purg_var_df.index.duplicated()]


## interpret variant tuples —— i.e. determine corresponding F/A genotypes at each mismatch site ##

# make a list of the changes that are present
all_changes = sorted(
						set(
								change 
										for entry in purg_var_df.index.levels[0] 
											for changes in eval(entry) 
											 	for change in changes
							)
)

change_count = len(all_changes)
if change_count != EXPECTED_SITE_COUNT:
	print(
			f"WARNING: Found {change_count} rather than the expected {EXPECTED_SITE_COUNT}", 
			"polymorphic sites!\a",
			file=sys.stderr
	)

# convert 'Tailed' level of column headings from Y(es) or N(o) into Tailed or tailess
#	try/except because this throws an error if no 'Tailed' level
try:
	purg_var_df.rename(
						axis='columns', 
						level='Tailed',
						mapper=lambda tail: 'Tailed' if tail.lower() in 'yes' else 'tailless',
						inplace=True
	)
except KeyError:
	pass

# change first (0th) level of index to a space separated string encoding whether each polymorphic 	# TODO: maybe add string version rather than convert —- ¿might make later processing easier?
#	site has the reference (F) allele or the templated alternate (A) allele
# purg_var_df.rename(
# 					axis='index', level=INDEX_NAMES[0], 
# 					mapper=lambda fset: " ".join(
# 													"A" if str(change) in fset else "F"  
# 														
# 													for change in all_changes
# 					),
# 					inplace=True
# )


# add level (to right/inside of 'Variant') named 'Variant string' parsed to F|A at each position #
index_df = purg_var_df.index.to_frame()

index_df.insert(	
				loc=1,
				column=f"{INDEX_NAMES[0]} string",
				value=[
						" ".join(
									"A" if str(change) in fset else "F"  
										
									for change in all_changes
						)
						
						for fset in purg_var_df.index.get_level_values(INDEX_NAMES[0])
				]
)

purg_var_df.index = pandas.MultiIndex.from_frame(index_df)

#### start MM2_fix stuff ####
# LOGIC: most likely a PCR/sequencing mutant if MM1 & MM3 agree and MM2 doesn't
# 	potential problem is that I'm not allowing for possibility that non-concordant 2nd position is correct 
#		 
#	TODO: ¿should be checking that frequency is higher for the 'correct' version before merging? 
MM2merged_purg_var_df = purg_var_df.copy()

MM2merged_purg_var_df.index = pandas.MultiIndex.from_frame(
													purg_var_df
														.index
														.to_frame()
														.replace(
																	regex={
																			'^F A F':'F F F',
																			'^A F A':'A A A'
																	},
														)
)

MM2merged_purg_var_df = MM2merged_purg_var_df.groupby(level='Variant string', sort=False).sum()

# add back MultiIndex levels that were lost in groupby and put all this back into purg_var_df <UGLY<<<<<<<<<<< TODO: make this less weird ... somehow
# 	√ TODO: this destroys the MultiIndex, which fucks up later steps FIXED
purg_var_df = (MM2merged_purg_var_df
											.reindex(
														purg_var_df.index,
														axis='index',
														level='Variant string'
											)
											.dropna(how='all')
)
#### end MM2_fix stuff ####

# original version of above renaming — made a tuple rather than space separated string
# purg_var_df.rename(
# 						axis='index', level=0, 
# 						mapper=lambda fset: tuple(
# 													"A" if str(change) in fset else "F"  
# 																	for change in all_changes
# 						),
# 						inplace=True
# )

# save relabeled variant frequencies
write_df_to_tsv(purg_var_df, out_file_base, index_label=None)


# original = broken and inelegant
# 	orig_mm2_merge_df = purg_var_df.copy()
# 
# 	variant_indexer = purg_var_df.index.get_level_values('Variant string').str.contains
# 
# 	FtoA_index = purg_var_df.index[variant_indexer('^F A F')].get_level_values('Variant string')
# 	AtoF_index = purg_var_df.index[variant_indexer('^A F A')].get_level_values('Variant string')
# 
# 	dubious_index = FtoA_index.append(AtoF_index)
# 	# good_index = (FtoA_index
# 	# 						.str.replace(
# 	# 										'F A F',
# 	# 										'F F F',
# 	# 										n=1,
# 	# 										regex=False
# 	# 						)
# 	# 						.append(
# 	# 									AtoF_index
# 	# 										.str.replace(
# 	# 														'A F A',
# 	# 														'A A A',
# 	# 														n=1,
# 	# 														regex=False
# 	# 										)
# 	# 						)
# 	# )
# 	good_index = (dubious_index
# 								.str.replace(
# 												'F A F',
# 												'F F F', 	# {'F A F':'F F F', 'A F A':'A A A'},
# 												n=1,
# 												regex=False
# 								)
# 								.str.replace(
# 												'A F A',
# 												'A A A',
# 												n=1,
# 												regex=False
# 								)
# 
# 	)
# 
# 	dubious_to_good = dict( zip(dubious_index, good_index) )
# 
# 	orig_mm2_merge_df.index = pandas.MultiIndex.from_frame(
# 															orig_mm2_merge_df
# 																.index
# 																.to_frame()
# 																.replace(dubious_to_good)
# 	)
# 
# 	orig_mm2_merge_df = orig_mm2_merge_df.groupby(level='Variant string', sort=False).sum()
# 
# 	# LOGIC: most likely a PCR/sequencing mutant if MM1 & MM3 agree and MM2 doesn't
# 	# 	potential problem is that I'm not allowing for possibility that non-concordant 2nd position is correct 
# 	#		 
# 	#	TODO: ¿should be checking that frequency is higher for the 'correct' version before merging? 
# 	mm2_merge_df = purg_var_df.copy()
# 
# 
# 	mm2_merge_df.index = pandas.MultiIndex.from_frame(
# 														mm2_merge_df
# 															.index
# 															.to_frame()
# 															.replace(
# 																		regex={
# 																				'^F A F':'F F F',
# 																				'^A F A':'A A A'
# 																		},
# 															)
# 	)
# 
# 	mm2_merge_df = mm2_merge_df.groupby(level='Variant string', sort=False).sum()
# 
# 	write_df_to_tsv(mm2_merge_df, f"mm2_merged_{out_file_base}", index_label=None)

# 	FtoA_seqs = mm2_merge_df.index[mm2_merge_df.index.get_level_values('Variant string').str.contains('^F A F')].get_level_values('Variant string')
# 	AtoF_seqs = mm2_merge_df.index[mm2_merge_df.index.get_level_values('Variant string').str.contains('^A F A')].get_level_values('Variant string')
# 	dubious_seqs = FtoA_seqs.append(AtoF_seqs)
# 	good_seqs = (FtoA_seqs
# 					.str.replace(
# 									'F A F',
# 									'F F F',
# 									n=1,
# 									regex=False
# 					)
# 					.append(
# 								AtoF_seqs
# 									.str.replace(
# 													'A F A',
# 													'A A A',
# 													n=1,
# 													regex=False
# 									)
# 					)
# 	)
# 
# 	dubious_to_good = dict( zip(dubious_seqs, good_seqs) )
# 
# 	# this doesn't work — ¿indicies don't match?
# 	mm2_merge_df[mm2_merge_df.index.get_level_values('Variant string').isin(good_seqs)] += mm2_merge_df[mm2_merge_df.index.get_level_values('Variant string').isin(dubious_seqs)]
# 
# 	dubious_df = mm2_merge_df[mm2_merge_df.index.get_level_values('Variant string').isin(dubious_seqs)]
# 	good_df = mm2_merge_df[mm2_merge_df.index.get_level_values('Variant string').isin(good_seqs)]
# 
# 	dubious_df_index = dubious_df.index.to_frame()
# 	dubious_df_index['Putative variant'] = dubious_df_index['Variant string'].transform(dubious_to_good)
# 
# 	dubious_df.index = pandas.MultiIndex.from_frame(dubious_df_index)
# 
# 	#dubious_df_index.replace(dubious_to_good)
# 	#dubious_df.index = pandas.MultiIndex.from_frame(dubious_df_index)
# 	dubious_df.index = pandas.MultiIndex.from_frame(dubious_df_index.replace(dubious_to_good))
# 		dubious_to_good = dict( zip( dubious_seqs.tolist(), good_seqs.tolist() ) )
#		dubious_to_good2 = dict(
# 								zip(
# 										dubious_seqs, 
# 										good_seqs 
# 								) 
# 		)
# 
# 		mm2_merge_df = mm2_merge_df[
# 		mm2_merge_df[mm2_merge_df.index.get_level_values('Variant string').contains()]
# 
# 					[
# 						for 
# 					]
# 		dubious_sequences = mm2_merge_df.index[mm2_merge_df.index.get_level_values('Variant string').str.contains('^A F A|^F A F')].get_level_values('Variant string').tolist()
# 
# 		for bad_seq in dubious_sequences:
# 	
# 			good_seq = bad_seq.replace(
# 	
# 			mm2_merge_df.xs(	
# 
# 		mm2_merge_df.index[mm2_merge_df.index.get_level_values('Variant string').str.contains('^A F A|^F A F')].tostring()
# 		mm2_merge_df.index.get_level_values('Variant string').str.contains('^A F A|^F A F')
# 		mm2_merge_df.index.get_level_values('Variant string').tolist()


#### suboptimal way to do this ... ###


# this is used below to make collapsed_df
#	goal is to add frequencies for sequences that are probably mutant derivatives to their parents
# row_corr_df = purg_var_df.T.corr()
# 
# nonidentical_row_corr_df = row_corr_df.where(row_corr_df < 1)
# 
# good_corrs_df = row_corr_df.where((row_corr_df < 1) & (row_corr_df > CORR_CUTOFF))
# 
# labeled_min_df = best_indices_to_header(purg_var_df, row_corr_df)	#, CORR_CUTOFF=0.5, KEEP_CUTOFF=0.13)
# 
# write_df_to_tsv(labeled_min_df, f"labeled_consolidated_{out_file_base}", index_label=None)

# old version, was leaving "best2" as tuples, which was difficult to interpret in spreadsheet & 
#	somehow messed up gt summary below ...
# OLD_head_df['best2'] = [
# 					tuple(
# 							sorted(
# 									(
# 										tup[0] 
# 												for tup in purg_var_df
# 															.nlargest(2, col, keep='all')
# 															.index.to_list()
# 												if purg_var_df.loc[tup, col].squeeze() > KEEP_CUTOFF				# weirdness going on here: more rational versions df.loc[tup, col] & df[col][tup] don't work! —— now trying .squeeze()
# 									),
# 									reverse=True
# 							)
# 					) 
# 					for col in purg_var_df.columns
# ]


# label each column with a guess as to what the genotype probably is
#	LOGIC: pick two most abundant sequences and only keep those with > KEEP_CUTOFF
#	PROBLEMS: 	1) what is the "correct/best" value for KEEP_CUTOFF
#				2) handle ambiguous cases – a) no sequence above KEEP_CUTOFF
#											b) >2 sequences pass (i.e. high score(s) with >1 seq.)
# TODO: functionalize?
# TODO: somehow deal with ambiguous results? —— probably best to consult before trying to do this
# TODO: determine if this is really working
# previous version, doesn't work with modified index (added stringified level)
# head_df = purg_var_df.columns.to_frame()
# head_df['best2'] = [
# 					"   ".join(
# 								sorted([			# simplifies results & loses relative frequency
# 											tup[0] 
# 										
# 											for tup in purg_var_df
# 														.nlargest(2, col, keep='all')
# 														.index.to_list()
# 										
# 											if purg_var_df.loc[tup, col].squeeze() > KEEP_CUTOFF				# weirdness going on here: more rational versions df.loc[tup, col] & df[col][tup] don't work! —— now trying .squeeze()
# 										],
# 										reverse=True
# 								)
# 					)
# 					for col in purg_var_df.columns
# ]
columns_df = purg_var_df.columns.to_frame()

# add (what will become the) 'best2' level
columns_df['best2'] = [
						"   ".join(			# sort simplifies results & loses relative frequency
									sorted([
												tup[1] 
										
												for tup in purg_var_df
															.nlargest(2, col, keep='all')
															.index
															.to_list()
										
												if purg_var_df.loc[tup, col].squeeze() > KEEP_CUTOFF	# weirdness going on here: more rational versions df.loc[tup, col] & df[col][tup] don't work! —— now trying .squeeze() √
											],
											reverse=True
									)
						)
						for col in purg_var_df.columns
]

# works, but just gives A-proportions as strings
# def duplex_inferred_from_freqs(df, col, change, F_cutoff=KEEP_CUTOFF, A_cutoff=1-KEEP_CUTOFF):
# 	"""
# 	"""
# 	return str(
# 													purg_var_df[
# 																	purg_var_df
# 																		.index
# 																		.get_level_values('Variant')
# 																		.str
# 																		.contains(
# 																					str(change), 
# 																					regex=False
# 																		) 
# 													][col].sum()
# 	)

def A_proportion(df, col, change):
	"""
	"""
	return df[
				df.index
					.get_level_values('Variant')
					.str
					.contains(
								str(change), 
								regex=False
					) 
	][col].sum()
	
#TOLERANCE = 0.15

columns_df['Sangerish'] = [
							repr(
								pandas.cut(
											[
												A_proportion(purg_var_df, col, change)
												for change in all_changes
											],
											bins=[	#bins, # or just use 5 (makes 5 equal sized bins)
													0.0,
													TOLERANCE,
													0.5 - TOLERANCE,
													0.5 + TOLERANCE,
													1.0 - TOLERANCE,
													1.0
											],
											labels=GENOTYPES,
											include_lowest=True
								)#.to_string()
							)
							for col in purg_var_df.columns
]

# columns_df['Sangerish'] = [
# 								"\t".join(
# 											str(
# 													purg_var_df[
# 																	purg_var_df
# 																		.index
# 																		.get_level_values('Variant')
# 																		.str
# 																		.contains(
# 																					str(change), 
# 																					regex=False
# 																		) 
# 													][col].sum()
# 											)
# 											for change in all_changes
# 								)
# 							for col in purg_var_df.columns
# ]


# 		bins = [0, F_CUTOFF, FA_LOW_CUTOFF, FA_HIGH_CUTOFF, A_CUTOFF, 1.0]
# 
# 		genotype_df = pandas.concat([ 
# 										pandas.cut(
# 													combined_purged_df[column], 
# 													bins=bins, 
# 													labels=GENOTYPES,
# 													include_lowest=True
# 										)
# 										for column in combined_purged_df.columns
# 									],
# 									axis=1
# 		)
# 
# 
# 		## check if there is any ambiguous data ##
# 		# ALTERNATIVE: if not genotype_df.isin(['F','F|A','A']).all().all():
# 		if genotype_df.isin( ["F|A?", "F?|A"] ).any().any():
# 			print("\nWARNING: Some sites were ambiguous — cutoffs may need adjustment.\n\a", file=sys.stderr)
# 
# 
# 		print(	
# 				f"{SPACER_HEAVY}\n"
# 				f"Categorical statistics for '{KEEP_ROW}' data by sequence position:\n"
# 				f"{SPACER_LIGHT}\n"
# 				f"{genotype_df.describe(include='all')}\n"
# 				f"{SPACER_HEAVY}\n"
# 				f"Sequence patterns for '{KEEP_ROW}' data:\n"
# 				f"{SPACER_LIGHT}\n"
# 				f"{genotype_df.apply(tuple, axis=1).value_counts()}\n"
# 				f"{SPACER_HEAVY}\n"
# 		)
# 
# 
# 		## save categorical version of data ##
# 		write_df_to_tsv(genotype_df, out_file_base + "_genotypes.tsv")
# 
# 
# 		## get frequencies of colony genotypes ##
# 		genotype_freq_df = pandas.concat([
# 											genotype_df[column].value_counts()
# 											for column in genotype_df.columns
# 										],
# 										axis=1,
# 										sort=True
# 		)

from itertools import product
def sum_crossprod_nlargest(df, n, col, keep='all'):
	return sum([a*b for a,b in product(df[col].nlargest(n, keep=keep),repeat=2)])

# TODO: handle "ties"
def scores_v0(df, col, n2):
	bestn2_prop = df.nlargest(n2, col, keep='first')[col]
	ultimate_prob = bestn2_prop[0] ** 2
	penultimate_prob = bestn2_prop[0] * bestn2_prop[1]
	
	return ultimate_prob, penultimate_prob, sum([a*b for a,b in product(bestn2_prop,repeat=2)])

def scores(df, col, n2, strand_sep=STRAND_SEP):
	bestn2_prop = df[col].droplevel(['Variant','Number of switches']).nlargest(n2, keep='first')
	
	proportion_ratio = bestn2_prop[0] / bestn2_prop[1]
	
	ultimate_prob = bestn2_prop[0] ** 2
	penultimate_prob = bestn2_prop[0] * bestn2_prop[1] #* 2	# reasonable to multiply by 2?
	
	top2_prop_sum = sum(bestn2_prop[0:2])
	
	probability_ratio = ultimate_prob / penultimate_prob
	
	top2_prob_total = sum_crossprod_nlargest(df, 2, col)
	 
	
	# prob vs. next best choice
	ultimate_vs_top2_prob_total = ultimate_prob / top2_prob_total
	penultimate_vs_top2_prob_total = 2 * penultimate_prob / top2_prob_total		# ¿reasonable to multiply by 2?
	
	if ultimate_vs_top2_prob_total > 3 * penultimate_vs_top2_prob_total:		# ¿too lax? (using: 3 *)	# ¿too strict?
		
		best = f"{bestn2_prop.index[0]}"
		abs_prob = ultimate_prob
		rel_prob = ultimate_vs_top2_prob_total									# simpler, but similar idea: ultimate_prob / (ultimate_prob + penultimate_prob)
		
		support = bestn2_prop[0]
		
	else:
		
		best = f"{bestn2_prop.index[0]}{strand_sep}{bestn2_prop.index[1]}"
		abs_prob = penultimate_prob
		rel_prob = penultimate_vs_top2_prob_total								# simpler, but similar idea: penultimate_prob / (ultimate_prob + penultimate_prob)
		
		support = top2_prop_sum
		
	return best, abs_prob, rel_prob, support, proportion_ratio, probability_ratio	#ultimate_prob, penultimate_prob, sum([a*b for a,b in product(bestn2_prop,repeat=2)])

# def scores_vector(df, col, n2=9999):
# 	sorted_col = df[col].nlargest(n2)
# 	ultimate_prob = bestn2_prop[0] ** 2
# 	penultimate_prob = bestn2_prop[0] * bestn2_prop[1] #* 2	# reasonable to multiply by 2?
# 	
# 	top2_prob_total = sum_crossprod_nlargest(df, 2, col)
# 	 
# 	
# 	# prob vs. next best choice
# 	ultimate_vs_top2_prob_total = ultimate_prob / top2_prob_total
# 	penultimate_vs_top2_prob_total = 2 * penultimate_prob / top2_prob_total
# 	
# 	if ultimate_vs_top2_prob_total > penultimate_vs_top2_prob_total:
# 		
# 		abs_prob = ultimate_prob
# 		rel_prob = ultimate_prob / top2_prob_total
# 		support = bestn2_prop[0]
# 		
# 	else:
# 		
# 		abs_prob = penultimate_prob
# 		rel_prob = 2 * penultimate_prob / top2_prob_total
# 		support = sum(bestn2_prop[0:2])
# 		
# 	return ultimate_prob, penultimate_prob, sum([a*b for a,b in product(bestn2_prop,repeat=2)])



#	limit {n->1/0} sum([a*b for a,b in product(purg_var_df.nlargest(n, purg_var_df.columns[0], keep='all')[purg_var_df.columns[0]],repeat=2)]) = 1
#	sum([a*b for a,b in product(purg_var_df.nlargest(2, purg_var_df.columns[0], keep='all')[purg_var_df.columns[0]],repeat=2)])
# score for the 'best2' 				TODO: finish this ...
def column_score(n=3, columns_df=columns_df, score_fn=scores):
	"""
	"""
	columns_df[f"score_{n}"] = [
									(
										score_fn(
													purg_var_df,
													col,
													n
										)
									)
									for col in purg_var_df.columns
	]
	
	return columns_df

columns_df = column_score(columns_df=columns_df, n=3)

# 
# columns_df['score_3'] = [
# 						(
# 							scores(
# 									purg_var_df,
# 									col,
# 									3
# 							)
# # 												purg_var_df
# # 															.nlargest(10, col, keep='all')[col]
# # 										
# # 												
# 						)
# 						for col in purg_var_df.columns
# ]
# 
# 
# columns_df['score_4'] = [
# 						(
# 							scores(
# 									purg_var_df,
# 									col,
# 									4
# 							)
# 						)
# 						for col in purg_var_df.columns
# ]
# 
# columns_df['score_5'] = [
# 						(
# 							scores(
# 									purg_var_df,
# 									col,
# 									5
# 							)
# 						)
# 						for col in purg_var_df.columns
# ]


purg_var_df.columns = pandas.MultiIndex.from_frame(columns_df)

write_df_to_tsv(purg_var_df, f"labeled_{CUTOFF_LABEL}_{out_file_base}", index_label=None)




# summarize genotypes found for each genotype tail combination #
# TODO: just write to a file rather than saving in a dictionary?
 
# # originally separated 	
# genotype_summaries_dict = {
# 							f"{lib}_{gt}_{'Tailed' if tail.lower() in 'yes' else 'tailless'}"			#f"{gt}_{'Tailed' if tl=='Yes' else 'tailless'}":			# key
# 							:
# 							df.columns.get_level_values('best2').value_counts()		# value (Series)
# 							
# 							for (lib, gt, tail), df in purg_var_df.groupby(
# 																				axis='columns', 
# 																				level=[
# 																						'Library',
# 																						'Genotype',
# 																						'Tailed'
# 																				]
# 							) 
# }
# 
# genotype_summary_df = pandas.DataFrame(genotype_summaries_dict).dropna(how='all')

# TODO: get column headings and put them in desired order 

# grouping_dict = {
# 					level: sorted(purg_var_df.columns.unique(level))
# 					
# 					for level in LEVELS
# }
# grouping = [
# 				(lib, gt, tail)
# 				
# 				for (lib, gt, tail), df 
# 					in purg_var_df.groupby( axis='columns', level=LEVELS )
# ]
# 
# from operator import itemgetter
# grouping.sort(key=itemgetter(2,1,0))
grouped_purg_var_df = purg_var_df.groupby(axis='columns', level=LEVELS)

# generate an ordering of tuples corresponding to values of LEVELS (levels used in reverse order)
#	use to reorder columns —— ordering based on LEVELS(in reverse order)
# grouping = sorted(
# 					(
# 						(lib, gt, 'Tailed' if tail.lower() in 'yes' else 'tailless')
# 				
# 						for (lib, gt, tail), df in grouped_purg_var_df
# 					),
# 					key=itemgetter(*reversed(range(len(LEVELS))))
# )

# get an ordered version of the column headings for levels used in groupby
grouping = sorted( 
					grouped_purg_var_df.indices.keys(), 
					key=itemgetter(*reversed(range(len(LEVELS))))
)

# 
# # manual reordering for clarity 																	#NB: probably needs to be deleted later <<<<<<<<<<<<<<<<<< !!!!!!!!!
# grouping[2:] = reversed(grouping[2:])

# TODO: Extract per position information 															<<<<<<<< in progress <<<<<<<<
#	may need to happen before/without converting tuples to space separated strings of [F|A]
#[purg_var_df.columns.get_level_values(INDEX_NAMES[0]).str.contains(change) for change in all_changes]	# list comprehension generating Boolean 

#purg_var_df.loc[ [purg_var_df.columns.get_level_values(INDEX_NAMES[0]).str.contains(change) for change in all_changes] ]

# actual frequency of A allele at each templated position (assumes all sequences are real)
# intended for comparison with the same type of data calculated for the 'best 2' (genotype_summary_df)
#	(see: per_position_Afreq_from_GTs_df)
per_position_Afreq_df = pandas.DataFrame(
										{
											(pos+1, change)
											:
											purg_var_df.loc[ 
																purg_var_df
																	.index
																	.get_level_values('Variant')
																	.str
																	.contains(
																				str(change), 
																				regex=False
																	) 
															].sum(axis='index')

											for pos, change in enumerate(all_changes)
										}
).T

per_position_Afreq_df.index.names = ['Mismatch #', 'Change tuple']

# describe fails, but mean works ... WTF?
#output = per_position_freq_df.groupby(level=LEVELS, axis='columns').describe()

mean_per_position_Afreq_df = per_position_Afreq_df.groupby(level=LEVELS, axis='columns').mean()
mean_per_position_Afreq_df = mean_per_position_Afreq_df.reindex(columns=grouping)

sd_per_position_Afreq_df = per_position_Afreq_df.groupby(level=LEVELS, axis='columns').std()
sd_per_position_Afreq_df = sd_per_position_Afreq_df.reindex(columns=grouping)


mean_per_position_Afreq_df.plot(
									kind='bar', 
									title='Proportion (mean ± sd) of A at each position', 
									figsize=(17,11), 
									rot=0,
									yerr=sd_per_position_Afreq_df, 
									xlim=(-1,7),
									ylim=(0,1),
)

save_show_plot(f"mean_per_position_Afreq_{out_file_base}.png")

###### ¿¿¿¿TODO: redo above, but utilizing the best/top2 choices???? ######

# mean_per_position_Afreq_df.T.hist(
# #									title='Proportion (mean ± sd) of A at each position', 
# #									by=mean_per_position_Afreq_df.columns,
# 									figsize=(17,11), 
# 									xrot=0,
# 									sharex=True,
# 									sharey=True,
# 									bins=20,
# #									yerr=sd_per_position_Afreq_df, 
# #									xlim=(-1,7)
# )
# 
# mean_per_position_Afreq_df.plot(
# 									kind='bar', 
# 									title='Proportion (mean ± sd) of A at each position', 
# 									figsize=(17,11), 
# 									rot=0,
# #									subplot=True,
# 									sharex=True,
# 									sharey=True,
# 									yerr=sd_per_position_Afreq_df, 
# 									xlim=(-1,7)
# )

# # counts occurences of inferred genotypes 
# # TODO: somehow deal with ambiguous results? —— probably best to consult before trying to do this
# genotype_summary_df = pandas.DataFrame(
# 										{
# 											(
# 												lib, 
# 												gt,
# 												'Tailed' if tail.lower() in 'yes' else 'tailless'
# 											)
# 											:
# 											df.columns
# 												.get_level_values('best2')
# 												.value_counts()
# 											
# 											for (lib, gt, tail), df in grouped_purg_var_df
# 										},
# 										# order columns
# 										columns=pandas.MultiIndex.from_tuples(grouping)
# )

genotype_summary_df = pandas.DataFrame(
										{
											levels
											:
											df.columns
												.get_level_values('best2')
												.value_counts()
											
											for levels, df in grouped_purg_var_df
										},
										# order columns
										columns=pandas.MultiIndex.from_tuples(grouping)
)



# name column levels
genotype_summary_df.columns.names = LEVELS
#genotype_summary_df.rename_axis(columns={'Variant': 'Combined'}, inplace=True)

# reorder rows —— puts rows with larger values at top
genotype_summary_df = genotype_summary_df.reindex(
													index=genotype_summary_df
															.max(1)
															.sort_values(ascending=False)
															.index,
)
# 
# # alternative .index[0][0] -> .index.get_level_values('Variant')
# # 	NB: fails for heteroduplexes, must be a better way
# get_changes = lambda seq: purg_var_df[
# 										purg_var_df
# 											.index
# 											.get_level_values('Variant string') == seq 
# ].index[0][0]
# 
# genotype_summary_df.index = pandas.MultiIndex.from_tuples([(get_changes(seq), seq) for seq in genotype_summary_df.index])


# def old_decompose(seq_string):
# 	"""
# 	Separate a 'best2' sequence into components.
# 	
# 	
# 	"""
# 	seqs_as_tuples = [ 
# 						seq.split()
# 						
# 						for seq in seq_string.split('   ')
# 	]
# 		
# 	
# 	return [
# 				frozenset(chars)
# 				for chars in zip(*seqs_as_tuples)
# 	]	

# record 'haplotype' frequencies
# TODO: decide if this should be done after adding adding index levels breaking down 'haplotypes' to separated 'SNP' alleles
write_df_to_tsv(
					genotype_summary_df,
					f"genotype_summary_{CUTOFF_LABEL}_{out_file_base}",
					index_label='Haplotype(s)'#None
)


 
# ### ¿temp? stuff for comparing only some libraries ###
# to_use = ['FAtailWT_S3', 'msh6FAtail_S2', 'msh2FAtail_S4']
# 
# my_barh=genotype_summary_df[to_use].plot.barh(figsize=(17,13), xticks=[*range(20)])

# in progress plotting stuff #
my_barh=genotype_summary_df.plot.barh(
										figsize=(17,13),
										color=[],
# 										xticks=range(int(genotype_summary_df.max().max())),#, xticks=[*range(20)]) 						# TODO: either make xticks fit input data rather than fixed value or allow default?
# 										rot=90,
)
plt.suptitle(f"Horizontal bar plot of estimated duplex frequencies where {CUTOFF_LABEL}")									# TODO: make figsize a default (or better yet depend on screen size ...)
plt.subplots_adjust(bottom=0.05, left=0.20, top=0.95, right=0.98)
plt.gca().invert_yaxis()
save_show_plot(f"Duplex_frequencies_horizontal_bar_{CUTOFF_LABEL}_{out_file_base}.png")

###NEW###
my_barh=genotype_summary_df.loc[genotype_summary_df.max(1) > 9].plot.barh(							# use .div(0.96) to plot as a percentage rather than a count
																			figsize=(17,13),
																			fontsize=14,
# 																			color=[
# 																					'dodgerblue',
# 																					'lime',
# 																					'magenta',
# 																					'blue',
# 																					'green',
# 																					'crimson'
# 																			],
# 																			xticks=range(int(genotype_summary_df.max().max())),#, xticks=[*range(20)]) 						# TODO: either make xticks fit input data rather than fixed value or allow default?
# 																			rot=90,
)
plt.suptitle(f"Horizontal bar plot of estimated duplex frequencies where {CUTOFF_LABEL}")									# TODO: make figsize a default (or better yet depend on screen size ...)
plt.subplots_adjust(bottom=0.05, left=0.23, top=0.95, right=0.98)
plt.gca().invert_yaxis()
save_show_plot(f"Top_duplex_frequencies_horizontal_bar_{CUTOFF_LABEL}_{out_file_base}.png")

# add tuple of tuples where each tuple corresponds to whether that position is being called  
#	as F, A, or F|A (heteroduplex) -> separate into individual levels of index
# should be able to use this to calculate the per position frequencies ...
genotype_summary_df.index = pandas.MultiIndex.from_tuples(
															[
																tuple(decompose(seq) + (seq,)) 
																for seq in genotype_summary_df
																			.index
															]
)

# ¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿and then what???????????????????????????
#values_df = pandas.DataFrame(columns=genotype_summary_df.columns)

# sets order of 'allele' occurrence — this is to match that used by Elena
alleles = (	
			frozenset('F'),
			frozenset('A'),
			frozenset(('F','A')),
)

# genotype/allele at each position — NB: based on assignments done using KEEP_CUTOFF
genotype_by_position_df = pandas.DataFrame(
											{
												( pos+1, '|'.join(sorted(chars, reverse=True)) )
												: 
												genotype_summary_df.loc[
																		genotype_summary_df
																			.index
																			.get_level_values(pos)
																			== 
																			chars
												].sum()
												
												for pos in range(change_count)
													for chars in (alleles)
											},
											dtype=int,
).T

genotype_by_position_df.index.names = ['Mismatch #', 'Allele']

# record 'allele' frequencies at each position
# TODO: decide if this should be done after adding adding index levels breaking down 'haplotypes' to separated 'SNP' alleles
write_df_to_tsv(
					dataframe=genotype_by_position_df,
					out_file_name=f"genotype_by_position_{CUTOFF_LABEL}_{out_file_base}",
					index_label=None
)



# 
# #genotype_by_position_df = genotype_by_position_df.unstack().reindex(columns=["|".join(allele) for allele in alleles], level='Allele')
# 
# (genotype_by_position_df
#     .unstack()
#     .reindex(
#     			columns=["|".join(allele) for allele in alleles],
#     			level='Allele'
#     )
#     .groupby(
#                             level=LEVELS,
#                             axis='columns'
#     )
#     .plot(
# #    		by=LEVELS,
# 			title=f"Genotype frequencies by position for {LEVELS}",
#  			kind='bar',
# 			figsize=(17,11), 
# 			rot=0,
# 			color=['b', 'r', 'g'],
# 			stacked=True,
# 			legend=False,
# # 			colormap='jet_r',
# #			subplot=True,
# 			sharex=True,
# 			sharey=True,
# 	)		
# )
# #     .bar(
# #                     figsize=(17,11), 
# #                     rot=0,
# #                     color=['b', 'r', 'g'],
# #                     stacked=True,
# #                     legend=False,
# # #					colormap='jet_r',
# # #					subplot=True,
# #                     sharex=True,
# #                     sharey=True,
# #     )
# )
# 
# # genotype_by_position_df.groupby(by='Mismatch #').plot(kind='bar', stacked=True, sharey=True, figsize=(17,11))
# # genotype_by_position_df.groupby(by='Mismatch #').plot(kind='bar', stacked=True, sharey=True, subplots=True, figsize=(17,11))
# # 
# # genotype_by_position_df.unstack().reindex(columns=["|".join(allele) for allele in alleles], level='Allele').plot(kind='bar', stacked=True, sharey=True, by=LEVELS, subplots=True, figsize=(17,11))
# 
# 


## above done separately instead of all together ... ##


# 
# # reshape df (and reorder) then groupby on LEVELS
# genotype_by_position_grouper = (
# 								genotype_by_position_df
# 									.unstack()
# 									.reindex(
# 												columns=["|".join(allele) for allele in alleles],
# 												level='Allele'
# 									)
# 									.groupby(
# 												level=LEVELS,
# 												axis='columns'
# 									)
# )
# 
# # plot & save bar chart for each subset of data
# for (levels, df) in genotype_by_position_grouper:
# 	
# 	strain_name = '_'.join(levels)
# 	
# 	df.plot(
# 				title=f"Genotype frequencies by position for {strain_name} where {CUTOFF_LABEL}",
# 				kind='bar',
# 				figsize=(17,11), 
# 				rot=0,
# 				color=['b', 'r', 'g'],
# 				stacked=True,
# 				legend=False,
# 	)
# 	
# 	save_show_plot(f"Genotype_frequencies_by_position_for_{strain_name}_{CUTOFF_LABEL}_{out_file_base}.png")


# redo above but with grouper done separately from reshaping ...
# reshape df (and reorder)
genotype_by_position_df = (
							genotype_by_position_df
									.unstack()
									.reindex(
												columns=[
															"|".join(
																		sorted(chars, reverse=True)
															) 
															for chars in alleles
												],
												level='Allele',
									)
)

# groupby LEVELS
genotype_by_position_grouper = genotype_by_position_df.groupby(
																level=LEVELS,
																axis='columns',
)


# TODO: repeat this using all data rather vs. best2 data (for Jim) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< *******!!!!!!
# plot & save bar chart for each subset of data
for (levels, df) in genotype_by_position_grouper:
	
	strain_name = '_'.join(levels)
	
	df.plot(
				title=f"Genotype frequencies by position for {strain_name} where {CUTOFF_LABEL}",
				kind='bar',
				figsize=(17,11), 
				rot=0,
				color=['b', 'r', 'g'],																# TODO: make into a default rather than hiding here ... or better yet use a dictionary to make order of colors correspond to genotype order?
				stacked=True,
				legend=False,
	)
	
	save_show_plot(f"Genotype_freq_by_position_for_{strain_name}_{CUTOFF_LABEL}_{out_file_base}.png")

# def xs_df(key, df=genotype_by_position_df, axis='columns', level='Allele'):							# <<<<<<<< now improved & (almost) always aborts ...) mysteriously this throws a KeyError the first time the script is run, but then works if rerun (happens everytime a new version is saved ...)
# 	"""
# 	Wrapper for xs applied to genotype_by_position_df.
# 	
# 	Input: key — value of level 'Allele'.
# 	
# 	Output: "slice" of df for axis & level where value is key.
# 	
# 	Requires: pandas.
# 	"""
# 	return df.xs(
# 					key=key,
# 					axis=axis,
# 					level=level,
# 	) 
	# xs_gt_by_pos_at_Allele = lambda key: genotype_by_position_Afreq_df.xs(
	# 																		key=key,
	# 																		axis='columns',
	# 																		level='Allele',
	# ) 

# df showing average frequency of 'A' at each position in each subset of data
#	for comparison with mean_per_position_Afreq_df
# TODO: automate comparisons?


df=genotype_by_position_df
axis='columns'
level='Allele'

genotype_by_position_Afreq_df = (
									( 
											2 * df.xs('A', axis=axis, level=level)
											+ df.xs('F|A', axis=axis, level=level)					# mysteriously this throws a KeyError the first time the script is run, but then works if rerun (happens everytime a new version is saved ...) —— update: actually works for KEEP_CUTOFF in [0.16, 0.18] ## and now back to not working (except for 0.16 ... ARGHHHH
									)
									/ ( 2 * genotype_by_position_grouper.sum() )									
)

write_df_to_tsv(
					genotype_by_position_Afreq_df, 
					f"mean_per_position_Afreq_from_gt_freqs_for_{CUTOFF_LABEL}_{out_file_base}", 
					index_label=None
)


genotype_by_position_Afreq_df.plot(
									kind='bar', 
									title=f"Proportion of A at each position where {CUTOFF_LABEL}", 
									figsize=(17,11), 
									rot=0,
#									yerr=sd_per_position_Afreq_df, 
									xlim=(-1,7),
									ylim=(0,1),
)

save_show_plot(f"mean_per_position_Afreq_from_gt_freqs_for_{CUTOFF_LABEL}_{out_file_base}.png")


diff_df = mean_per_position_Afreq_df - genotype_by_position_Afreq_df
diff_vector = diff_df.stack([0,1,2])																#TODO: Fix this — ¿assumes three level MultiIndex for columns?

def count_haplotypes(duplex_seqs):
	"""
	"""
	try:
		result = len(duplex_seqs)//13
	except:
		if duplex_seqs is numpy.NaN:
			# expected result if no sequence at all
			result = 0
		else:
			# shouldn't ever happen so raise hell (nonsense result)
			result = -666
	
	return result
		

#
duplex_haplotype_vector = genotype_summary_df.index.get_level_values(-1).map(count_haplotypes)



def mean_genotype_freq(
						keys,
						levels=('Genotype', 'Tailed', 'Allele'),
						axis='columns',
						df=genotype_by_position_df
):
	"""
	"""
	return (
				df.xs(
						keys,
						axis=axis,
						level=levels
				)
				.mean()
				/
				df.xs(
						keys[:2],
						axis=axis,
						level=levels[:2]
				)
				.sum(axis=axis)
				.mean()
	).squeeze()

#msh2_tailless_freq_AF = mean_genotype_freq(('msh2', 'tailless', 'F|A'))
#msh6_tailless_freq_AF = mean_genotype_freq(('msh6', 'tailless', 'F|A'))


print(
		"Quantification of distance between measured and predicted (from assigned genotypes)"
		f" 'A' frequencies for {KEEP_CUTOFF}:\n"
		f"\tSUM of DIFFERENCES: {diff_df.sum().sum()}\n"
		f"\tSUM of ABSOLUTE DIFFERENCES: {diff_df.abs().sum().sum()}\n"
		f"\tSUM of SQUARED DIFFERENCES: {diff_df.pow(2).sum().sum()}\n"
		f"\tMEAN ± SD of DIFFERENCES: {diff_vector.mean()} ± {diff_vector.std()}\n"
		f"\tFIVE NUMBER SUMMARY of DIFFERENCES: {diff_vector.quantile([0,0.25,0.5,0.75, 1.0])}\n"
		f"\tMAD of DIFFERENCES: {diff_vector.mad()}\n"
		f"\tSKEW of DIFFERENCES: {diff_vector.skew()}\n"
		f"\tKURTOSIS of DIFFERENCES: {diff_vector.kurt()}\n"
		f"\tNumbers of distinct strands per duplex [min, .. , max]: {sorted(duplex_haplotype_vector.unique())}\n"
#		f"\tMean FA frequency for msh2 tailess: {mean_genotype_freq(('msh2', 'tailless', 'F|A'))}\n"
#		f"\tMean FA frequency for msh6 tailess: {mean_genotype_freq(('msh6', 'tailless', 'F|A'))}\n"
		f"\tFrequencies of top 10 duplexes:\n{genotype_summary_df.droplevel([0,1,2,3,4,5,6]).loc[['F F F F F F F', 'F F F F F F F   A A A A A A A', 'A A A A A A A', 'F F F A A A A   A A A A A A A', 'F F F A A A A', 'F F F F F F F   F F F A A A A', 'F F F A A A A   F F A A A A A', 'F F A A A A A   A A A A A A A', 'F F A A A A A', 'F F F F F F F   F F A A A A A']]}\n"
		# TODO: plot these (as proportions) ... see also scoring only version ...
)
# 
# # temp testing
# genotype_summary_df.droplevel([0,1,2,3,4,5,6]).loc[['F F F F F F F', 'F F F F F F F   A A A A A A A', 'A A A A A A A', 'F F F A A A A   A A A A A A A', 'F F F A A A A', 'F F F F F F F   F F F A A A A', 'F F F A A A A   F F A A A A A', 'F F A A A A A   A A A A A A A', 'F F A A A A A', 'F F F F F F F   F F A A A A A']]


# 
# per_position_Afreq_from_GTs_df.index.names = ['Mismatch #', 'Change tuple']



# for (levels, df) in (
# 						genotype_by_position_df
# 							.unstack()
# 							.reindex(
# 										columns=["|".join(allele) for allele in alleles],
# 										level='Allele'
# 							)
# 							.groupby(
# 													level=LEVELS,
# 													axis='columns'
# 							)
# ):
# 	df.plot(
# 	#    		by=LEVELS,
# 				title=f"Genotype frequencies by position for {levels}",
# 				kind='bar',
# 				figsize=(17,11), 
# 				rot=0,
# 				color=['b', 'r', 'g'],
# 				stacked=True,
# 				legend=False,
# 	# 			colormap='jet_r',
# 	#			subplot=True,
# # 				sharex=True,
# # 				sharey=True,
# 	
# 	)
# 	
# 	save_show_plot(f"Genotype_frequencies_by_position_for_{levels}_{CUTOFF_LABEL}_{out_file_base}.png")



# genotype_by_position_df.unstack().sort_index(
# 												axis='columns',	# was: 0
# 												level=-1,	# was: None
# 												ascending=False,	# was: True
# 												inplace=False,
# 												kind='quicksort',
# 												na_position='last',
# 												sort_remaining=False,	# was: True
# 												by=None
# )
												

#genotype_by_position_df.unstack().groupby(level=LEVELS, axis='columns').plot.bar(stacked=True, sharey=True) #, subplots=True


# #Available color maps:
# 				0			1		2			3		4		5			6		7		8		9			10			11		12			13		14		15			x16		17			18			19		20		21			22			23			?24		25		26			27			28			29			30
# COLORMAPS = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']
# 
# grpby = genotype_by_position_df.unstack().groupby(level=LEVELS, axis='columns')
# 
# def color_test():
# 	color_num = 0
# 	while color_num < len(COLORMAPS):
# 		for (lib, gt, tail), df in grpby:
# 			colors = COLORMAPS[color_num]
# 			color_num += 1
# 			print(f"using: {colors}")
# 			df.plot.bar(
# 							figsize=(17,11), 
# 							rot=0,
# 							colormap=colors,
# 							stacked=True
# 			)
# 			plt.show()

	


#genotype_by_position_df.plot.bar(subplots=True, sharey=True)

# genotype_by_position_df = pandas.DataFrame(
# 											index={
# 												(pos+1, chars)
# 												: 
# 												genotype_summary_df.loc[
# 																			genotype_summary_df
# 																				.index
# 																				.get_level_values(pos)
# 																				== 
# 																				chars
# 												].sum()
# 												
# 												for pos in range(change_count)
# 													for chars in (frozenset('F'),frozenset(('F','A')),frozenset('A'))
# 											},
# 											columns=genotype_summary_df.columns
# )
# 
# genotype_by_position_df = pandas.DataFrame(
# 											index={
# 												(pos+1, chars)
# 												: 
# 												genotype_summary_df.loc[
# 																			chars in 
# 																			genotype_summary_df
# 																				.index
# 																				.get_level_values(pos)
# 												].sum()
# 												
# 												for pos in range(change_count)
# 													for chars in (('F',),('F','A'),('A'))
# 											},
# 											columns=genotype_summary_df.columns
# )
# 
# 								{
# 									(pos, chars): 
# 									
# 									genotype_summary_df.loc[indx]
# 									for indx in genotype_summary_df.index
# 										for pos, chars in enumerate(decompose(indx))
# 								},
# 								columns=genotype_summary_df.columns
# )
# values_df = pandas.DataFrame(
# 								{
# 									(pos, chars): 
# 									
# 									genotype_summary_df.loc[indx]
# 									for indx in genotype_summary_df.index
# 										for pos, chars in enumerate(decompose(indx))
# 								},
#								columns=genotype_summary_df.columns
# )
# for indx in genotype_summary_df.index:
# 	for pos, chars in enumerate(decompose(indx)):
# 		values_df.loc[(pos,chars)] += genotype_summary_df.loc[indx]




# # counts occurences of inferred genotypes 
# # TODO: somehow deal with ambiguous results? —— probably best to consult before trying to do this
# genotype_summary_df = pandas.DataFrame(
# 										{
# 											(
# 												lib, 
# 												gt,
# 												'Tailed' if tail.lower() in 'yes' else 'tailless'
# 											)
# 											:
# 											df.columns
# 												.get_level_values('best2')
# 												.value_counts()
# 											
# 											for (lib, gt, tail), df 
# 												in purg_var_df.groupby(
# 																		axis='columns', 
# 																		level=LEVELS
# 												)
# 											
# 										},
# 										#columns=[]
# 										
# )
# 
# 
# # generate an ordering of tuples corresponding to values of LEVELS (levels used in reverse order)
# grouping = sorted(
# 					(
# 						(lib, gt, 'Tailed' if tail.lower() in 'yes' else 'tailless')
# 				
# 						for (lib, gt, tail), df 
# 							in purg_var_df.groupby(axis='columns', level=LEVELS)
# 					),
# 					key=itemgetter(*reversed(range(len(LEVELS))))
# )
# 
# # reorder columns —— ordering based on LEVELS(in reverse order)
# genotype_summary_df = genotype_summary_df.reindex(
# 													columns=pandas.MultiIndex.from_tuples(grouping),
# 													#copy=False
# )
# 
# # reorder columns; name column levels; reorder rows #
# genotype_summary_df.columns.names = LEVELS
# #genotype_summary_df.rename_axis(columns={'Variant': 'Combined'}, inplace=True)
# 
# # reorder rows —— puts rows with larger values at top
# genotype_summary_df = genotype_summary_df.reindex(
# 													index=genotype_summary_df
# 															.max(1)
# 															.sort_values(ascending=False)
# 															.index,
# )
# 





# # think about this, but probably need to add something when creating ... ^^^^^^^^^^^^
# lvls_direction = [('Library', False), ('Genotype', True)]
# for lvl, ascend in lvls_direction:
# 	genotype_summary_df = genotype_summary_df.reindex(
# 														level=lvl,
# 														columns=genotype_summary_df
# 																	.columns
# 																	.unique(lvl)
# 																	.sort_values(ascending=ascend),
# 	)
# 
# lvl, ascend = 'Library', False; genotype_summary_df = genotype_summary_df.reindex( level=lvl, columns=genotype_summary_df.columns.unique(lvl).sort_values(ascending=ascend) )
# 
# # this doesn't quite seem to work
# genotype_summary_df = genotype_summary_df.reindex(
# 													level='Genotype',
# 													columns=genotype_summary_df
# 																.columns
# 																.unique('Genotype')
# 																.sort_values(ascending=False),
# )
# 
# # this is AFU —— also not clear how to presort when creating genotype_summary_df ...				<<<<<<<<< figure this out ....
# genotype_summary_df = genotype_summary_df.reindex(
# 													level='Genotype',
# 													columns=genotype_summary_df
# 																.columns
# 																.get_level_values('Genotype')
# 																.sort_values(ascending=False),
# )
# 
# 




# split up into subplots by column
#my_barh=genotype_summary_df.plot.barh(figsize=(17,14), subplots=True, sharey=True, fontsize=8); plt.suptitle("Horizontal bar plot of estimated duplex frequencies") ; plt.subplots_adjust(bottom=0.05, left=0.14, top=0.95, right=0.98) ; plt.gca().invert_yaxis()

# not split up
#my_barh=genotype_summary_df.plot.barh(figsize=(17,13), fontsize=12); plt.suptitle("Horizontal bar plot of estimated duplex frequencies") ; plt.subplots_adjust(bottom=0.05, left=0.2, top=0.95, right=0.98) ; plt.gca().invert_yaxis()
 
# # table experiments:
# ugly failure and includes figure ... (plt.axis('off') before my_table is same)
# my_table=genotype_summary_df.T.plot(figsize=(17,13), fontsize=12, table=True, axes=False); plt.axis('off'); plt.suptitle("Horizontal bar plot of estimated duplex frequencies") ; plt.subplots_adjust(bottom=0.5, left=0.06, top=1.0, right=0.98)
# # almost works, but I can't get the text big enough ...
# plt.axis('off'); my_table=plt.table(cellText=genotype_summary_df.values, rowLabels=genotype_summary_df.index, colLabels=genotype_summary_df.columns, fontsize=60); plt.suptitle("Horizontal bar plot of estimated duplex frequencies") ; plt.subplots_adjust(bottom=0.5, left=0.06, top=1.0, right=0.98)


## get top sequences —— currently all those with at least 1 occurence of a proportion > 0.1 ##
# minimum proportion to count
# PROPORT_CUTOFF = 0.05	# was 0.1
# # minimum counts to include
# COUNT_CUTOFF = 1

# for each row count the number of columns that contain a proportion greater than 0.1 
freq_seq_count = purg_var_df[purg_var_df >= PROPORT_CUTOFF].count(axis='columns')
# get the 'sequences' (names) that have counts (see previous line) greater than 0 
top_seqs = freq_seq_count[freq_seq_count >= COUNT_CUTOFF].index.tolist()

#top10_variants_df = purg_var_df[purg_var_df.rank(ascending=False) <= 10].dropna(how='all')


# 1st: get genotype prevalence by colony data
# 2nd: get genotype prevalence by genotype/tailing combination
#	for each genotype/tailing combination
#		summarize results for (all observed | all observed with freq > x% | all observed for gt&tail with freq/rank > x? | top n? ) sequence patterns
#		Y: heatmap

#	for group in purg_var_df.groupby(axis='columns', level=[1,2]): group[0], group[1].agg([numpy.mean, numpy.std], axis=1)
#	for group in purg_var_df.groupby(axis='columns', level=[1,2]): if group[0][0]: print(f"{group[0]}\n{group[1].quantile(numpy.divide(range(5), 4), axis=1)}")
#plt.rcParams['font'] = "Courier"

 
# # ugly ...
# box_plots = purg_var_df.loc[top_seqs].groupby(axis='columns', level=['Genotype','Tailed']).quantile(QUARTILES, axis=1).plot(kind='box', figsize=(11,8.5), vert=False); plt.show()
 
# # FAILS ... doesn't seem to deal with subplots ... may be a version issue?
# box_plots = purg_var_df.loc[top_seqs].groupby(axis='columns', level=['Genotype','Tailed']).quantile(QUARTILES, axis=1).plot(kind='box', figsize=(11,8.5), rot=90, subplots=True); plt.show()
 
# # too crowded
# area_plots = purg_var_df.loc[top_seqs].groupby(axis='columns', level=['Genotype','Tailed']).plot(kind='area', figsize=(11,8.5), subplots=True, rot=90, legend=None, stacked=True, sharey=True, xticks=range(len(top_seqs))); plt.show()


# for (gt, tail), df in purg_var_df.groupby(axis='columns', level=['Genotype','Tailed']): 
for levels, df in purg_var_df.loc[top_seqs].droplevel([0,2]).groupby(
																		axis='columns', 
																		level=LEVELS
				# 														[
				# 																'Library',
				# 																'Genotype',
				# 																'Tailed'
				# 														]
): 
	
	strain_name = '_'.join(levels)			#"{lib}_{gt}_{'Tailed' if tail.lower() in 'yes' else 'tailless'}"		#f"{gt}_{'Tailed' if tail=='Yes' else 'tailless'}"
	
	# get rid of extraneous 'Number of switches' information (from vast)
#	df.index = df.index.droplevel(1)
 
# 	# narrow down to sequence patterns that: 
# 	#	make up > 5% of sequences in one sample with current genotype/tail combination
# 	#	occur with a median proportion > 1% of samples for current genotype/tail combination
#  	quartiles_df = df[ (df.max(1) > 0.05) & (df.median(1) > 0.01) ].quantile(QUARTILES, axis=1)
	quartiles_df = df.quantile(QUARTILES, axis=1)
		
#	temp_df = grp[1][(grp[1].max(1) > 0.05) & (grp[1].median(1) > 0.01)].quantile(numpy.divide(range(5), 4), axis=1)
# 	temp_df.rename(	
# 						axis='columns', level=0,
# 						mapper=lambda x: " ".join(x),
# 						inplace=True
# 		)
		
	quartiles_name = f"{strain_name}__quartiles_{CUTOFF_LABEL}_{out_file_base}"
	
	write_df_to_tsv(
						dataframe=quartiles_df, 
						out_file_name=quartiles_name, 
						index_label='Quartiles'
					)
	# boxplot (using quartiles) of data from current df
#	box_plot = temp_df.boxplot(rot=90, figsize=(11,8.5), showmeans=True)#, notch=True, bootstrap=1000) #, sym="", bottom=0.2, top=0.95)
#	plt.subplots_adjust(bottom=0.2, top=0.95)
	box_plot = quartiles_df.boxplot(figsize=(11,8.5), vert=False, showmeans=True)#, subplots=True, notch=True, bootstrap=1000) #, sym="", bottom=0.2, top=0.95)
	plt.subplots_adjust(bottom=0.05, left=0.16, top=0.95, right=0.96) #plt.subplots_adjust(bottom=0.05, left=0.2, top=0.95)
#	plt.show()
	plt.suptitle(f"Boxplot of allele frequencies within {strain_name} where {CUTOFF_LABEL}")
	plt.gca().invert_yaxis()
	save_show_plot(quartiles_name + ".png")
 	
# 	# figure out 
# 	# maybe useful??
# 	unstacked_area_plot=df.plot.area(figsize=(11,8.5), rot=90, stacked=False, legend=None, xticks=range(len(df.index)))
# 	plt.subplots_adjust(bottom=0.2, left=0.08, top=0.95, right=0.96) #plt.subplots_adjust(bottom=0.2, top=0.95)
# 	plt.suptitle(f"Unstacked area plot of allele frequencies within {strain_name} where {CUTOFF_LABEL}")
# 	save_show_plot(f"{strain_name}__unstacked_area_{CUTOFF_LABEL}_{out_file_base}.png")
# 	
# 	# maybe useful??
# 	stacked_area=df.plot.area(figsize=(11,8.5), rot=90, stacked=True, legend=None, xticks=range(len(df.index)))
# 	plt.subplots_adjust(bottom=0.2, left=0.06, top=0.95, right=0.96) #plt.subplots_adjust(bottom=0.2, top=0.95)
# 	plt.suptitle(f"Stacked area plot of allele frequencies within {strain_name} where {CUTOFF_LABEL}")
# 	save_show_plot(f"{strain_name}__stacked_area_{CUTOFF_LABEL}_{out_file_base}.png")
	
	my_barh=df.plot.barh(figsize=(11,8.5), legend=None)
	plt.suptitle(f"Horizontal bar plot of allele frequencies within {strain_name} where {CUTOFF_LABEL}")
	plt.subplots_adjust(bottom=0.05, left=0.16, top=0.95, right=0.96) #plt.subplots_adjust(bottom=0.05, left=0.2, top=0.95)
	plt.gca().invert_yaxis()
	save_show_plot(f"{strain_name}__horizontal_bar_{CUTOFF_LABEL}_{out_file_base}.png")

# 	# below don't work on actual data ... only seeing data for first sequence????
# 	#	TODO figure out what is actually being plotted!
# 	myplot=df.plot.density(figsize=(11,8.5))
# 	myplot=df.plot.kde(figsize=(11,8.5))									# v. similar to above ...
#	myplot=df.plot.hist(figsize=(11,8.5))									# wtf?
# 	
# 	#row_corr_df = df.T.corr() --- decided to use row_corr for entire dataset not just this subset ... TODO: implement this
# 	
# 	collapsed_df = df.copy()
# 	CORR_CUTOFF = 0.5
# 	
# 	# analyze each column (colony) separately ––– current version is quite slow ... TODO: speed up?
# 	for col in df.columns:
# 		# look at each row (sequence variant) starting with most common 
# 		for row in collapsed_df[col].sort_values(ascending=False).index:
# 			
# 			# get top two largest (largest will be self) correlations as a pandas Series
# 			#two_largest_corrs = row_corr_df[row].nlargest(2, row)
# 			# above fails, but below works —— bug in Series implementation of nlargest for whatever version of pandas I'm currently using ...
# 			two_largest_corrs = row_corr_df.nlargest(2, row)[row]
# 			
# 			# is this correlation better than the (arbitrary) cutoff value?
# 			if two_largest_corrs[1] > CORR_CUTOFF:
# 				# most correlated other row (sequence variant), putative PCR induced variant of row
# 				corr_row = two_largest_corrs.index[1]
# 				
# 				# is the putative PCR variant at least as complex as its supposed parent?
# 				if corr_row in collapsed_df.index:
# 					# check to make sure adding less frequent to more frequent sequence
# 					if collapsed_df.loc[row, col] < collapsed_df.loc[corr_row, col]:
# 						# switch row and corr_row so that sequence with larger proportion is being kept
# 						row, corr_row = corr_row, row
# 					# check that putative PCR variant at least as complex as its supposed parent
# 					if corr_row[1] >= row[1]: 	
# 						# add value of variant to putative parent
# 						collapsed_df.loc[row, col] += collapsed_df.loc[corr_row, col]	#[0]
# 						# replace variant proportion with NaN
# 						collapsed_df.loc[corr_row, col] = numpy.NaN
# 		
# 		#
# 		
# 	
# 	# purge all rows that now entirely NaN
# 	collapsed_df.dropna(how='all', inplace=True)
# 	
# 	# this still needs work, but it is getting there ... converts MultiIndex to Index, better to modify columns MultiIndex ...
# 	collapsed_df.append( pandas.DataFrame([[collapsed_df.nlargest(2, col)[col].index.to_flat_index() for col in df.columns]], columns=collapsed_df.columns, index=('best2',))  )
# 	# this just adds a new row
# 	collapsed_df.loc['best2'] = [collapsed_df.nlargest(2, col)[col].index.to_flat_index() for col in df.columns]
# 	collapsed_df.loc['best2'] = [list(collapsed_df.nlargest(2, col)[col].index) for col in df.columns]
# 	# gets rid of unneeded "number of switches" info
# 	collapsed_df.loc['best2'] = [[tup[0] for tup in list(collapsed_df.nlargest(2, col)[col].index)] for col in df.columns]
# 	
# 	#	pandas.DataFrame([[msh2_untailed.nlargest(2, col)[col].index.to_flat_index() for col in msh2_untailed.columns]], columns=msh2_untailed.columns,index=('best2',)) 
# ### eliminate weird duplicate sequences (errors from vast) ###





column_corr_df = purg_var_df.corr()
#column_corr.droplevel(1).droplevel(1, axis=1)

# use cutoffs to segment and then groupby corr?


# statistics for inter and intra group correlations ...
# identify weird samples?

# 
# row_corr = purg_var_df.T.corr()
# row_corr = row_corr.droplevel(1).droplevel(1, axis=1)
# n=5	# NB: based on data this should get 
# for col in row_corr.columns: 
# 	print(f"'{col}' {n}-largest correlations\n======================================\n{row_corr.nlargest(n, col)[col]}\n")
# 
# 
# top_row_corr = purg_var_df.T[top_seqs].corr()
# n=5
# for col in top_row_corr.columns: 
# 	print(f"'{col}' {n}-largest correlations\n======================================\n{top_row_corr.nlargest(n, col)[col]}\n")
# 
#row_corr = row_corr.reindex(row_corr.rank().sort_values(by=[('F F F F F F F', '0'),('A A A A A A A', 1)]).index, axis='index')
#row_corr.reindex(row_corr.sort_values(ascending=False).index, axis='index')

#"F F F F F F F"
#"A A A A A A A"


### pick best allele(s) to allow calling genotype associated with each colony
# 
# # current rules for addressing weirdness at MM2 (probably due to Taq polymerase messiness):
# #	• corr > 0.5 between variant and parent
# #	• if variant position 'conflicts' with L & R positions then non-conflicting sequence = parent
# #	• if L & R disagree then use most common allele with a ? (cf using max & count > 0.1)
# variant_parent_dict = {
# 						'F A F F F F F':'F F F F F F F', 	# corr > 0.7
# 						'F A F A A A A':'F F F A A A A',	# corr > 0.9
# 						'F A F F F F A':'F F F F F F A',	# corr > 0.8
# 						'A F A A A A A':'A A A A A A A',	# corr > 0.5
# 						'F A F F F A A':'F F F F F A A',	# corr > 0.6
# 						'F A F A A F A':'F F F A A F A',	# corr > 0.8
# # 						'F F A A A A A':'F F? A A A A A',	# corr > 0.7 with 'F A A A A A A'
# # 						'F A A A A A A':'F F? A A A A A',	# corr > 0.7 with 'F F A A A A A'
# # 						'A A F F F F F':'A A? F F F F F',	# corr > 0.6 with 'A F F F F F F'
# # 						'A F F F F F F':'A A? F F F F F'	# corr > 0.6 with 'A A F F F F F'
# }
# 
# simp_var_df = purg_var_df.copy()
# 
# # this works, but seems kludgy
# #	¿better to do this from top to bottom in each column so that don't get wrong direction conversions?
# #		would need to redo variant_parent_dict to use that most effectively ... automatable?
# #	alternatives: do this on groups sorting in between additions
# #				  do this on each colony separately, working upward
# #				  when add rows make new row with origin row names combined ...
# for variant, parent in variant_parent_dict.items():
# 	simp_var_df.loc[parent] = simp_var_df.loc[parent].values + simp_var_df.loc[variant].values		#### this wont work the way I've made variant_parent_dict!!! ####
# 	simp_var_df.drop(variant, level=0, inplace=True)
# 
# # **** how about step by step merge starting with lowest and working until only two sequences left???????
# # analyze each colony separately – determine actual state of starting DNA duplex in SSA region
# top2_variants_df = simp_var_df[simp_var_df.rank(ascending=False) <= 2].dropna(how='all')
# 
# 
# # retry ...
# variant_parent_dict3 = {
# 						('F A F F F F F', 1):('F F F F F F F', 0), 	# corr > 0.7
# 						('F A F A A A A', 2):('F F F A A A A', 1),	# corr > 0.9
# 						('F A F F F F A', 2):('F F F F F F A', 1),	# corr > 0.8
# 						('A F A A A A A', 2):('A A A A A A A', 1),	# corr > 0.5
# 						('F A F F F A A', 2):('F F F F F A A', 1),	# corr > 0.6
# 						('F A F A A F A', 3):('F F F A A F A', 2),	# corr > 0.8
# # 						('F F A A A A A', 1):('F F? A A A A A', 1),	# corr > 0.7 with 'F A A A A A A'
#  						('F A A A A A A', 1):('F F A A A A A', 1),	# corr > 0.7 with 'F F A A A A A'
# # 						('A A F F F F F', 1):('A A? F F F F F', 1),	# corr > 0.6 with 'A F F F F F F'
# # 						('A F F F F F F', 1):('A A? F F F F F', 1)	# corr > 0.6 with 'A A F F F F F'
# }
# 
# 
# simp_var_df3 = purg_var_df.loc[top_seqs].copy()
# for variant, parent in variant_parent_dict.items():
#     simp_var_df3.loc[parent] = simp_var_df3.loc[parent].values + simp_var_df3.loc[variant].values#### this wont work the way I've made variant_parent_dict!!! ####
#     simp_var_df3.drop(variant, inplace=True)





#### pseudo pseudocode ####
# Goal is to make a reasonable guess as to what the GT was in the parental cell of the colony #
#	¿address ties by sorting by overall proportion (for particular gentoype/tail dataset)?
#	for each of the (top n?) values in column starting with largest
#		add values from less common possible variants
#	score / compare top 2 values — ratio | difference | ??
#	if exceed some threshold then just use sequence associated with top value
#	otherwise use top two
# 	associate score with this call for later use?



# for column in top10_variants_df.columns:
# 	print( top10_variants_df[column][top10_variants_df[column] > 0].sort_values(ascending=False) )
	
# purg_var_df[purg_var_df > 0].groupby(axis='columns', level=[1,2])


# purg_var_df.apply(	
# 							sort,
# 							axis='index'
# )
	


# 			change_length = sum(
# 									change[2] if change[1]=='D' else len(change[2]) 
# 									for change in changes
# 							)



##				all_changes = sorted( set( change for entry in purg_var_df.index for changes in eval(entry[0]) for change in changes ) )
# all_changes = sorted(
# 						set(	
# 								change 
# 										for entry in purg_var_df.index 
# 											for changes in eval(entry[0]) 
# 												for change in changes
# 							)
# 					)

# alternatives ...
##				all_changes2 = sorted( set( change for entry in purg_var_df.index.levels[0] for changes in eval(entry) for change in changes ) )
# slightly faster version ...

# purg_var_df.index.levels[0]
# 
# 
# purg_var_df.insert( loc=2, column='Alleles', value=purg_var_df.index.levels[0] )
# 

#for change in all_changes:
	




##				all_changes2b = sorted( set( changes for entry in purg_var_df.index.levels[0] for changes in eval(entry) ) )
# this gives all present combinations ...
# all_changes2b = sorted(
# 						set(
# 								changes 
# 										for entry in purg_var_df.index.levels[0] 
# 											for changes in eval(entry) 
# 							)
# 					)

# all_changes3 = sorted(
# 						set(
# 								*changes 
# 										for entry in purg_var_df.index.levels[0] 
# 											for changes in eval(entry) 
# 							)
# 					)

### below is doing something weird ...
#sorted(	set(	change for entry in purg_var_df.index.levels[0] for change in zip(*eval(entry))	)	)


## for each 'frozenset' step through and create an F or A label for each known site ... ##
new_level = purg_var_df.index.levels[0].map( lambda fset: ["A" if str(change) in fset else "F"  for change in all_changes] )


# ## output statistics for columns containing data ##
# print(	
# 		f"{SPACER_HEAVY}\n"
# 		f"Statistics for '{KEEP_ROW}' data by sequence position:\n"
# 		f"{SPACER_LIGHT}\n"
# 		f"{combined_purged_df.describe()}\n"
# 		f"{SPACER_HEAVY}\n"
# )
# 
# ## save results for columns containing data ##
# write_df_to_tsv(combined_purged_df, out_file_base + "_frequencies.tsv")
# 
# 
# ## report if there isn't data for the expected number of sequence positions ##
# found_count = len(combined_purged_df.columns)
# if found_count != EXPECTED_SITE_COUNT:
# 	print(
# 			f"WARNING: Found {found_count} rather than the expected {EXPECTED_SITE_COUNT}"
# 			"polymorphic sites!\a",
# 			file=sys.stderr
# 	)
# 
# 
# ## make, save, and display histogram of frequencies ##
# freq_plot = combined_purged_df.hist(bins=20, figsize=(11,8.5), sharey=True, range=(0,1))
# plt.suptitle('A-allele frequency at each variable position')
# save_show_plot(out_file_base + "_frequencies.png")
# 
# ## convert colony frequencies into genotypes ##
# # bin boundaries made using arbitrary CUTOFF values 
# bins = [0, F_CUTOFF, FA_LOW_CUTOFF, FA_HIGH_CUTOFF, A_CUTOFF, 1.0]
# 
# genotype_df = pandas.concat([ 
# 								pandas.cut(
# 											combined_purged_df[column], 
# 											bins=bins, 
# 											labels=GENOTYPES
# 								)
# 								for column in combined_purged_df.columns
# 							],
# 							axis=1
# )
# 
# 
# ## check if there is any ambiguous data ##
# # ALTERNATIVE: if not genotype_df.isin(['F','F|A','A']).all().all():
# if genotype_df.isin( ["F|A?", "F?|A"] ).any().any():
# 	print("\nWARNING: Some sites were ambiguous — cutoffs may need adjustment.\n\a", file=sys.stderr)
# 
# 
# print(	
# 		f"{SPACER_HEAVY}\n"
# 		f"Categorical statistics for '{KEEP_ROW}' data by sequence position:\n"
# 		f"{SPACER_LIGHT}\n"
# 		f"{genotype_df.describe(include='all')}\n"
# 		f"{SPACER_HEAVY}\n"
# 		f"Sequence patterns for '{KEEP_ROW}' data:\n"
# 		f"{SPACER_LIGHT}\n"
# 		f"{genotype_df.apply(tuple, axis=1).value_counts()}\n"
# 		f"{SPACER_HEAVY}\n"
# )
# 
# 
# ## save categorical version of data ##
# write_df_to_tsv(genotype_df, out_file_base + "_genotypes.tsv")
# 
# 
# ## get frequencies of colony genotypes ##
# genotype_freq_df = pandas.concat([
# 									genotype_df[column].value_counts()
# 									for column in genotype_df.columns
# 								],
# 								axis=1,
# 								sort=True
# )
# 
# ## save frequencies of colony genotypes ##
# write_df_to_tsv(genotype_freq_df, out_file_base + "_genotype_frequencies.tsv", index_label="Genotype")
# 
# 
# ### create plot showing different heteroduplex patterns ############################################
# # NB: problem with below approach is that it requires assuming that F will always be linked with F
# #		is that known to be true?
# # 	to obviate this, need to use data from variant_frequency report ...
# 
# ##genotype_df["combined"] = genotype_df.apply(lambda x: x.values, axis=1) # alternatively put this into a Series
# #genotype_df["combined"] = genotypes_df.apply(tuple, axis=1)
# 
# #genotype_representation_df = pandas.DataFrame(columns=["List of SNPs", "Frequency"])
# # store results
# duplex_freq_list = []
# 
# for genotype_tuple, count in genotype_df.apply(tuple, axis=1).value_counts().items():
# 	
# 	duplex = pandas.Series([[], []], index=["strand1","strand2"])
# 	# build representation of DNA duplex
# 	for element in genotype_tuple:
# 		
# 		# version 0
# 		chars = element.split('|')
# 
# 		try:
# 			duplex = duplex.add([ [chars[0]],[chars[1]] ])
# 		except IndexError:
# 			duplex = duplex.add([ [chars[0]],[chars[0]] ])
# 
# 	duplex_freq_list.append((duplex, count))
# 
# # needed?
# genotype_representation_df = pandas.DataFrame( duplex_freq_list, columns=["Duplex SNPs", "Frequency"])
# 
# ## save frequencies of colony genotypes ##
# write_df_to_tsv(genotype_representation_df, out_file_base + "_duplex_SNP_frequencies.tsv", index_label="Genotype")
# 
# # maybe what is useful is a format like this:	strand1:	F---F---F---F---F---F---F	69
# #												strand2:	F---A?--A---A---A---A---A
# # or maybe something that fits in a spreadsheet?
# #
# # TODO: better to apply/transform or something the elements in the genotype_representation_df
# 
# genotype_representation_df["Duplex"] =	[
# 											"\n".join(
# 														"---".join(strand) for strand in duplex
# 											) 
# 											for duplex, count in duplex_freq_list
# 										]
# 
# write_df_to_tsv(genotype_representation_df.loc[:,['Duplex', 'Frequency']], out_file_base + "_duplex_frequencies.tsv", index_label=None)
# 
# 
# 		
# # 	
# # 	genotype_representation_df["List of SNPs"] = 
# # 	genotype_representation_df["Frequency"] = 
# 		
# # 
# # 		# version 1
# # 		chars = element.split('|')
# # 
# # 		try:
# # 			duplex = duplex.add([ [char] for char in chars ])
# # 		except ValueError:
# # 			duplex = duplex.add([ [chars],[chars] ])
# # 		
# # 		# version 2
# # 		chars = element.split('|')
# # 
# # 		duplex.apply(lambda strand: strand.append(chars)
# # 
# # 
# # 		# version 3
# # 		chars = element.split('|')
# # 
# # 		if len(chars) == 2:
# # 			# heteroduplex site
# # 			duplex.
# # 		
# # 
# # 		# version 4
# # 		if "|" in element:
# # 			# heteroduplex site
# # 			strand1.append('F')
# # 			if "?" in element:
# # 				pass
# # 		
# # 		elif "F" in element:
# # 			# F (reference sequence) site
# # 			pass
# # 		
# # 		else:
# # 			# must be an A (cassette sequence) site
# # 			pass
# # 		
# 		
# 	