#!/usr/bin/python
# Scripts for converting .BED les to .interval_list les
# BED les are used by the coverage scripts while .interval_list are used by auto-classication

import argparse
import os, glob, shutil, subprocess
import warnings
#import fileinput
import datetime, time
import numpy
import sqlite3, csv, itertools
import math
from decimal import *
from suds.client import Client

from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from bs4 import BeautifulSoup
from lxml import etree
import xml.dom.minidom

from illuminate import InteropDataset
	
def main():
	parser = argparse.ArgumentParser(description='Python version of the legacy createReport.cygwin.sh')
	parser.add_argument('run_folder', help='Path to Run Folder')
	parser.add_argument('qscore', help = 'QScore threshold', type = int) 
	parser.add_argument('coverage', help = 'Coverage Threshold', type = int) # AW: Optional argument: --start option
	parser.add_argument('out_folder', help = 'Path to Analysis Out Folder For NGS-Analysis') 
	parser.add_argument('resource_folder', help = 'Path to folder containing pipeline scripts')
	parser.add_argument('--reanalysis', help = 'Indicate whether this is for reanalysis. (0:No, 1:Yes) [Default: 0]', choices=[0, 1], type=int, default = 0)
	parser.add_argument('--casefolder', help = 'Indicate case folder root to store/sort a copy of the NGS results. Leave empty if this option is not desired.')
	
	args = parser.parse_args()
	
	case_folder_suffix = ""
	
	if args.reanalysis == 1:
		print("Reanalysis mode.")
		case_folder_suffix = "/REANALYSIS_REPORTS"
	else:
		print("Analysis mode.")
		case_folder_suffix = "/[BREAST|COLON|UNSORTED_NGS_REPORTS]"
	
	if args.casefolder:
		print ("Copy of reports will be saved in " + args.casefolder + case_folder_suffix)
	
	# Set parameters
	warnings.simplefilter(action = "ignore", category = FutureWarning) # Disables Illuminate's sort()
	warnings.simplefilter(action = "ignore", category = RuntimeWarning) # Disables numpy's NaNMean
	
	run_id = os.path.basename(args.run_folder.rstrip("/"))
	sample_sheet = os.path.join(args.run_folder.rstrip("/"),"SampleSheet.csv")
	start_time = datetime.datetime.now()
	
	# CONSTANT - DBS / RESOURCES
	analysis_type_csv = os.path.join(args.resource_folder,"createReport_DBs","csv_db", "analysis_type.csv")
	metrics_db = os.path.join(args.resource_folder,"createReport_DBs","csv_db", "run_metrics.db")
	variants_db = os.path.join(args.resource_folder,"createReport_DBs","csv_db", "variants.db")	
	# CONSTANT - FOLDERS - DEFAULT BASED ON resource_folder, out_folder, and run_folder definitions
	script_folder = os.path.dirname(os.path.realpath(__file__))
	align_folder = os.path.join(args.run_folder.rstrip("/"),"Data","Intensities","BaseCalls","Alignment")
	settings_folder = os.path.join(args.out_folder, run_id, "Settings")
	analysis_settings_file = os.path.join(settings_folder, "analysis_settings.txt")
	report_folder = os.path.join(args.out_folder,run_id,"Reports")
	cnv_folder = os.path.join(align_folder, "cnv")
	seq_only_genelist = os.path.join(args.resource_folder,"createReport_DBs","analysis_type","seq_only.genelist") #NBN, BARD1, FANCC due to lack of MLPA
	cnv_only_genelist = os.path.join(args.resource_folder,"createReport_DBs","analysis_type","cnv_only.genelist") #Like EPCAM or GREM1 - Exclude from Variant Related Code

	# CONSTANT - SUFFIXES [Suffixes are file suffixes followed by SampleID but without the file extensions]
	bam_suffix = "-aln-pe-sorted-merged-dedup-recal"
	coverage_output_suffix = bam_suffix + "-covQ" + str(args.coverage)
	av_output_suffix = "_Annov_out.hg19_multianno"
	raw_vcf_suffix = "_raw_variants"
	
	# CONSTANT - DOCKER VOLUME MOUNTS - Need to mount root of: 1. Resource Folder, 2. Out Folder, 3. Run Folder, 4. Case Folder
	
	dockerVolume = ''
	docker_roots = []
	resource_path = args.resource_folder.lstrip(os.sep)
	resource_root = resource_path[:resource_path.index(os.sep)] if os.sep in resource_path else resource_path #basically finds first occurence of separater, and then do a 
	dockerVolume = '-v /' + resource_root + "/:/" + resource_root + '/'
	docker_roots.append(resource_root)
	
	out_path = args.out_folder.lstrip(os.sep)
	out_root = out_path[:out_path.index(os.sep)] if os.sep in out_path else out_path
	if out_root not in docker_roots:
		dockerVolume = dockerVolume + ' -v /' + out_root + "/:/" + out_root + '/'
		docker_roots.append(out_root)		
	
	run_path = args.run_folder.lstrip(os.sep)
	run_root = run_path[:run_path.index(os.sep)] if os.sep in run_path else run_path
	if run_root not in docker_roots:
		dockerVolume = dockerVolume + ' -v /' + run_root + "/:/" + run_root + '/'	
		docker_roots.append(run_root)
	
	if args.casefolder:
		case_path = args.casefolder.lstrip(os.sep)
		case_root = case_path[:case_path.index(os.sep)] if os.sep in case_path else case_path
		if case_root not in docker_roots:
			dockerVolume = dockerVolume + ' -v /' + case_root + "/:/" + case_root + '/'	
			docker_roots.append(case_root)
	
	"""
	dockerVolume='-v /media/:/media' # Most mounted media on VBox and the likes are in /media/
	
	For referemces only:
	dockerVolume='-v '+args.resource_folder+':'+args.resource_folder
	dockerVolume=dockerVolume+' '+'-v '+args.out_folder+':'+args.out_folder
	dockerVolume=dockerVolume+' '+'-v '+args.run_folder+':'+args.run_folder
	dockerVolume=dockerVolume+' '+'-v '+os.path.dirname(args.roi_bed_path)+':'+os.path.dirname(args.roi_bed_path)
	"""
	
	# CONSTANT - DOCKER COMMAND
	samtoolsCMD = "docker run --rm " + dockerVolume + " nderoo324/samtools"
	bedtoolsCMD = "docker run --rm " + dockerVolume + " nderoo324/bedtools bedtools"
	pythonCMD = "docker run --rm " + dockerVolume + " nderoo324/python"
	rscriptCMD = "docker run --rm " + dockerVolume + " nderoo324/r-base" #Rscript
	
	# GLOBAL VARIABLES - META_DICT & REPORT_QC_DICT
	software_version = "CreateReport V3"
	genome = "hg19"
	
	# INITIATING VARIABLES
	report_qc_dict = {}
	meta_dict = {}
	
	cnv_only_dict = get_gene_dict(cnv_only_genelist)
	seq_only_dict = get_gene_dict(seq_only_genelist)

	if not os.path.isdir(report_folder):
		os.makedirs(report_folder)
		os.makedirs(settings_folder)

	analysis_folder = os.path.join(align_folder,"Analysis")
	if not os.path.isdir(analysis_folder):
		os.makedirs(analysis_folder)
		os.makedirs(os.path.join(analysis_folder,"annovar"))
		os.makedirs(os.path.join(analysis_folder,"coverage"))
		
	run_date = run_id.split("_")[0]
	instrument_id = run_id.split("_")[1]
	flowcell_id = ""
	instrument = ''

	if instrument_id.startswith("M"):
		instrument = 'MiSeq'
		flowcell_id = run_id.split("-")[1]
	else: #instrument_id.startswith("N"):
		instrument = 'NextSeq'
		flowcell_id = run_id.split("_")[3]

	meta_dict['flowcell_id'] = flowcell_id
	meta_dict['run_date'] = run_date
	meta_dict['instrument_id'] = instrument_id
	meta_dict['report_date'] = time.strftime("%Y-%m-%d")
	meta_dict['coverage'] = args.coverage
	meta_dict['qscore'] = args.qscore
	meta_dict['genome'] = genome
	meta_dict['software_version'] = software_version
	meta_dict['variants_called'] = 0 # Defined based on filter_av in main()
	meta_dict['variants_reported'] = 0 # Defined based on variants extracted from filter_av in generate_report() with PASS flag
	meta_dict['variants_failed'] = 0 # To-Do: Need to define based on variants extracted from filter_av in generate_report() with non-PASS flag

	# =======
	# Start of script
	# =======
	 
	print("Checking required files...")
	requisite_flag = checkRequisites(args.run_folder, align_folder, instrument)
	
	if requisite_flag  == 0:
		exit()
		
	print("Configuring analysis from SampleSheet.csv...")
	
	# =======
	# Start recording QC in the background
	# =======
	if instrument == 'MiSeq':
		record_run_metrics (args.run_folder, align_folder, metrics_db)
	
	sample_dict = {}
	sample_dict = configure_analysis_by_sample_sheet(sample_sheet, analysis_type_csv, analysis_settings_file)

	print("Processing BAMS...")
	for sample_id in sample_dict:
		print ("\tFiltering by QScore and removing duplicates for " + sample_id + "...")
		filterBAM(align_folder, sample_id, args.qscore, coverage_output_suffix, bam_suffix, samtoolsCMD)
	
	print ("Calculating coverage and processing variants...")
	for sample_id in sample_dict:
		bed = sample_dict[sample_id]['bed']
		genelist = sample_dict[sample_id]['genelist']
		
		print ("\tCalculating coverage...")
		coverage_rscript = os.path.join(script_folder, "rscripts", "coverage_boxplot.R")
		[failed_exons_list, gene_exons_dict, exons_summary_csv] = get_coverage (align_folder, analysis_folder, sample_id, args.qscore, args.coverage, genelist, bed, coverage_rscript, coverage_output_suffix, bedtoolsCMD, rscriptCMD)
		sample_dict[sample_id]['failed_exons'] = failed_exons_list
		sample_dict[sample_id]['gene_exons'] = gene_exons_dict
		sample_dict[sample_id]['exons_summary_csv'] = exons_summary_csv
		
		print ("\tProcessing annotated variants...")
		av_info="DP,SB"
		av_geno="GT,GQ,AD"
		av_file = glob.glob(os.path.join(align_folder, sample_id + "*" + av_output_suffix + ".txt"))[0]
		flatten_av_file = os.path.join(analysis_folder, "annovar", sample_id + "_flatten_av.txt")
		flatten_av (av_file, flatten_av_file, av_info, av_geno) 
		# Extended Version (dbsnp and 1000g)
		#filter_av_field = "Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,avsnp144,1000g2015aug_all,ExAC_nontcga_ALL,Zygosity,QUAL,FILTER,DP,AD,GQ,GT,SB"
		filter_av_field = "Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,ExAC_nontcga_ALL,Zygosity,QUAL,FILTER,DP,AD,GQ,GT,SB"
		filter_av_file= os.path.join(analysis_folder, "annovar", sample_id + "_filter_av.txt")
		filter_av (flatten_av_file, filter_av_file, genelist, bed, filter_av_field, variants_db, cnv_only_dict) 

		print ("\tCalculating QC metrics for " + sample_id + "...")
		# Consolidate RunQC stats
		report_qc_dict[sample_id] = report_qc_miseq (args.run_folder, run_id, sample_id, instrument)
		
		# Update variants called for file - Note: The RAW VCF and the FLATTEN VCF now contains ALL the variants from panel. It is the FILTERED VCF that narrows it to the targeted region.
		with open(filter_av_file) as f:
			sample_dict[sample_id]['variants_called']  = sum(1 for _ in f) - 1 # The first line is header.

	for sample_id in sample_dict:
	
		print ("Generating Reports for " + sample_id + "...")
		coverage_folder = os.path.join(analysis_folder,"coverage")
		other_plots_folder = os.path.join(coverage_folder, "other_plots")
				
		if not os.path.isdir(other_plots_folder):
			os.makedirs(other_plots_folder)
	
		print ("\tGenerating Base Coverage plots for " + sample_id + "...")
		base_coverage_rscript = os.path.join(script_folder, "rscripts",  "base_vs_coverage_histo.R")
		base_coverage_log = os.path.join(analysis_folder,"coverage","other_plots",sample_id + "_base_vs_coverage_plot.log")
		command = " ".join([rscriptCMD, base_coverage_rscript, sample_id, analysis_settings_file , coverage_folder, other_plots_folder]) # requires sample_analysis_file and *all_only.hist
		subprocess.call(command, stdout=open(base_coverage_log, 'w'), stderr=subprocess.STDOUT, shell = True)
		
		print ("\tGenerating Alt Frequency plots for " + sample_id + "...")
		altfreq_coverage_rscript = os.path.join(script_folder, "rscripts", "altfreq_vs_coverage_plot.R")
		altfreq_coverage_log = os.path.join(analysis_folder,"coverage","other_plots",sample_id + "_altfreq_vs_coverage_plot.log")
		command = " ".join([rscriptCMD, altfreq_coverage_rscript, sample_id, analysis_settings_file , os.path.join(analysis_folder, "annovar"), other_plots_folder]) # requires sample_analysis_file and *varsummary.csv		
		subprocess.call(command, stdout=open(altfreq_coverage_log, 'w'), stderr=subprocess.STDOUT, shell = True)
	
		print ("\tExamining CNV Calls...")
		bed = sample_dict[sample_id]['bed']
		genelist = sample_dict[sample_id]['genelist']
	
		cnv_panel_csv = os.path.join(args.resource_folder, "bed_files", "cnv_panel_exons.csv")
		cnv_dict = parse_cnv(sample_id, cnv_folder, cnv_panel_csv, bed, genelist)
		cnv_plot_folder = os.path.join(cnv_folder, "exomedepth")
		cnv_plots_list = get_cnv_plots(sample_id, sample_dict[sample_id]['genelist'], cnv_plot_folder)
		
		# Create the PDF Reports
		print ("\tCreating PDFs for " + sample_id + "...")
		meta_dict['sample_id'] = sample_id
		meta_dict['manifest_id'] = sample_dict[sample_id]['manifest']
		meta_dict['experiment_id'] = sample_dict[sample_id]['experiment']
		meta_dict['investigator'] = sample_dict[sample_id]['investigator']
		meta_dict['variants_called'] = sample_dict[sample_id]['variants_called']
		av_file = os.path.join(analysis_folder, "annovar", sample_id + "_filter_av.txt")
		alt_plot = os.path.join(other_plots_folder, sample_id + "_altfreq_vs_coverage.png")
		base_plot = os.path.join(other_plots_folder, sample_id + "_bases_vs_coverage.png")
		boxplot_folder = os.path.join(coverage_folder, "gene_coverage", "boxplot")
		boxplots_list = glob.glob(boxplot_folder + "/" + sample_id + "*.png")
		#print sample_dict[sample_id]['exons_summary_csv']

		generate_report(args.resource_folder, script_folder, report_folder, args.casefolder, av_file, cnv_dict, report_qc_dict[sample_id], meta_dict, sample_dict[sample_id]['gene_exons'], sample_dict[sample_id]['failed_exons'], sample_dict[sample_id]['exons_summary_csv'], cnv_only_dict, seq_only_dict, alt_plot, base_plot, boxplots_list, cnv_plots_list, args.reanalysis)
		
		print ("\tReport for " + sample_id +  " created.")
	
	print "Script finished successfully"
	print "Script run time: " + str(datetime.datetime.now() - start_time)
			
######################################################################################################################################################
######################################################################################################################################################
########################################################################### END OF MAIN ##################################################################
######################################################################################################################################################
######################################################################################################################################################

def checkRequisites(run_folder, align_folder, instrument):

	requisite_found = 1
	
	if not os.path.isdir(align_folder):
		print ("Error: The required folder" + align_folder + "is not found.")
		print ("Please make sure that the desired alignment folder is named Alignment")
		requisite_found = 0
	
	if not os.path.isdir (os.path.join(run_folder,"InterOp")):
		print ("Error: The required folder 'InterOp' is not found in " + run_folder)
		requisite_found = 0
	
	if not os.path.exists (os.path.join(run_folder,"RunInfo.xml")):
		print ("Error: The required file 'RunInfo.xml' is not found in " + run_folder)
		requisite_found = 0
	
	#MiSeq Only
	if instrument == 'MiSeq':
		if not os.path.exists (os.path.join(align_folder,"Summary.xml")):
			print ("Error: The required le 'Summary.xml' is not found in " + align_folder)
			requisite_found = 0
	
	if not os.path.exists (os.path.join(run_folder,"SampleSheet.csv")):
		print ("Error: The required le 'SampleSheet.csv' is not found in " + run_folder)
		requisite_found = 0
	
	return requisite_found
	
def get_gene_dict(genelist):
	gene_dict = {}
	
	try:
		for row in open(genelist, 'rb'):
			row_fields = row.rstrip('\r\n').split('\t')
			gene_name = row_fields[0]
			accession = row_fields[1]
			gene_dict[gene_name] = accession
	except:
		pass

	return gene_dict
	
def configure_analysis_by_sample_sheet(sampleSheetName, analysisTypeCSV, analysisSettingsFile):
	sampleSheet = open(sampleSheetName, 'r')
	lines = sampleSheet.readlines()
	outFile = open(analysisSettingsFile, 'w')
	outFile.write("sample,analysis,caller,manifest,genelist,bed,experiment\n")
	
	manifestDict = {} # Holds the collection of manifests
	sampleSheetKeyDict = {} # Holds the column heading for sample sheet
	sampleDict = {} # Holds the detailed info for samples
	
	sampleSheetDict = {}
	analysisTypeDictReader = csv.DictReader(open(analysisTypeCSV))
	
	analysisTypeDict = []
	
	for analysis in analysisTypeDictReader:
		analysisTypeDict.append(analysis)
	
	sampleSheetDict['investigator'] = "-"
	sampleSheetDict['assay'] = "-"
	sampleSheetDict['description'] = "-"
	sampleSheetDict['experiment'] = "-"
	sampleSheetDict['indexi5'] = "-"
	sampleSheetDict['indexi7'] = "-"
	sampleSheetDict['caller'] = "GATK"
	sampleSheetDict['exclude'] = "-"
	sampleSheetDict['analysis_type'] = "-"
	
	i5present = True #sometimes it's missing
	
	for line in enumerate(lines):
		# Check if i5 is present. Changes the parsing of run qc data
		
		if "Investigator Name," in line[1]:
			sampleSheetDict['investigator'] = line[1].split(',')[1].rstrip('\r\n')
		elif "ExcludeRegion" in line[1]:
			if "," in line[1]:
				if line[1].rstrip('\r\n').split(',')[1]:
					sampleSheetDict['exclude'] =  "Error"
		elif "VariantCaller," in line[1]:
			sampleSheetDict['caller'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Assay Name," in line[1]:
			sampleSheetDict['assay'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Experiment Name," in line[1]:
			sampleSheetDict['experiment'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Description," in line[1]:
			sampleSheetDict['description'] = line[1].split(',')[1].rstrip('\r\n')
		elif "[Manifests]" in line[1]:
			i = 1
			# Assertion: At least ONE manifest is present, else there's an error in the sample sheet
			manifestLine = lines[line[0]+i]
			
			while manifestLine not in ['\n', '\r\n']:
				manifestID = manifestLine.split(',')[0]
				manifestFile = manifestLine.split(',')[1].rstrip('\r\n')
				manifestDict[manifestID] = {}
				manifestDict[manifestID]['manifest'] = manifestFile
				manifestDict[manifestID]['genelist'] = '-'
				manifestDict[manifestID]['bed'] = '-'
				manifestDict[manifestID]['analysis_type'] = '-'

				#print manifestID
				#print manifestFile
				
				# Associate analysis type with manifest
				for analysis in analysisTypeDict:
					#print analysis['manifest_name']
					if manifestFile == analysis['manifest_name']:
						#print manifestFile
						manifestDict[manifestID]['genelist'] = analysis['genelist']
						manifestDict[manifestID]['bed'] = analysis['bed']
						manifestDict[manifestID]['analysis_type'] = analysis['analysis_type']
						
				i=i+1
				manifestLine = lines[line[0]+i]
				# This line error checks for empty lines filled with commas
				manifestLine = manifestLine.lstrip(',')
	
		elif "Sample_ID" in line[1]:
			
			for idx,sampleKey in enumerate(line[1].rstrip('\r\n').split(',')):
				sampleSheetKeyDict[sampleKey] = idx
				
			if "I5_Index_ID" not in line[1]:
				i5present = False
					
			i = 1
			# Assertion: At least ONE manifest is present, else there's an error in the sample sheet
			sampleLine = lines[line[0]+i]
			
			while sampleLine not in ['\n', '\r\n']:
				
				sampleID = sampleLine.split(',')[0]
				
				if (i5present):
					sampleSheetDict['indexi5'] = sampleLine.split(',')[sampleSheetKeyDict['I5_Index_ID']] + ":" + sampleLine.split(',')[sampleSheetKeyDict['index2']]
				
				sampleSheetDict['indexi7'] = sampleLine.split(',')[sampleSheetKeyDict['I7_Index_ID']] + ":" + sampleLine.split(',')[sampleSheetKeyDict['index']]
	
				#NextSeq Manifest
				if len(manifestDict) == 0:
					manifestFile = sampleLine.split(',')[sampleSheetKeyDict['Description']].rstrip('\r\n')
	
                              	 	# Associate analysis type with manifest
                               		for analysis in analysisTypeDict:
                                        	if manifestFile == analysis['manifest_name']:
							sampleSheetDict['manifest'] = manifestFile
                                                	sampleSheetDict['genelist'] = analysis['genelist']
                                                	sampleSheetDict['bed'] = analysis['bed']
                                                	sampleSheetDict['analysis_type'] = analysis['analysis_type']

				#MiSeq Manifest
				else:
					manifestID = sampleLine.split(',')[sampleSheetKeyDict['Manifest']]
				
					if manifestID in manifestDict:
						sampleSheetDict['manifest'] = manifestDict[manifestID]['manifest']
						sampleSheetDict['analysis_type'] = manifestDict[manifestID]['analysis_type']
						sampleSheetDict['genelist']  = manifestDict[manifestID]['genelist']
						sampleSheetDict['bed']  = 	manifestDict[manifestID]['bed']
					else:
						sampleSheetDict['manifest'] = "-"
						sampleSheetDict['analysis_type'] = "-"
		
				sampleDict[sampleID] = sampleSheetDict.copy()			
				outFile.write(",".join([sampleID, sampleSheetDict['analysis_type'], sampleSheetDict['caller'],  sampleSheetDict['manifest'], sampleSheetDict['genelist'], sampleSheetDict['bed'],  sampleSheetDict['experiment']]) + "\n")

				i=i+1
				
				if line[0]+i == len(lines):
					break
					
				sampleLine = lines[line[0]+i]
				# This line error checks for empty lines filled with commas
				sampleLine = sampleLine.lstrip(',')

	outFile.close()
	
	return sampleDict

def record_run_metrics(run_folder, align_folder, metrics_db):
	
	runMetricArray = {}
	runID = os.path.basename(run_folder)
	fileMissing = 0

	if os.path.exists(os.path.join(align_folder, "aggregate.enrichment_summary.csv")):

		aggregateSummary = open(os.path.join(align_folder, "aggregate.enrichment_summary.csv"), "r")

		runMetricArray['ID'] = "'" + runID + "'"

		for row in aggregateSummary:
			rowArray = row.split(",")
			rowCat = rowArray[0]

			if rowCat == "Uniformity of coverage (Pct > 0.2*mean):":
				runMetricArray['UNIFORMITY'] = round(float(rowArray[1].rstrip("%")), 2)
			elif rowCat == "Percent duplicate paired reads:":
				runMetricArray['PCT_DUPLICATE'] = round(float(rowArray[1].rstrip("%")), 2)
			elif rowCat == "SNV Ts/Tv ratio:":
				runMetricArray['TI_TV_RATIO'] = round(float(rowArray[1]), 2)
			elif rowCat == "Total PF reads:":
				runMetricArray['TOTAL_PF_READS'] = int(rowArray[1])
			elif rowCat == "Percent aligned reads:":
				runMetricArray['PCT_ALIGNED'] = round(float(rowArray[1].rstrip("%")), 2)
			elif rowCat == "Fragment length median:":
				runMetricArray['MEDIAN_FRAGMENT_LENGTH'] = int(rowArray[1])
			elif rowCat == "Mean region coverage depth:":
				runMetricArray['MEAN_COVERAGE_DEPTH'] = round(float(rowArray[1].rstrip("%")), 2)
			elif rowCat == "Percent Q30:":
				runMetricArray['PCT_Q30'] = round(float(rowArray[1].rstrip("%")), 2)
			elif rowCat == "Read enrichment:":
				runMetricArray['PCT_ON_MANIFEST'] = round(float(rowArray[1].rstrip("%")), 2)
	else:
		fileMissing = 1

	if os.path.exists(os.path.join(align_folder, "EnrichmentStatistics.xml")):
		enrichStats = os.path.join(align_folder, "EnrichmentStatistics.xml")

		xmlDom = xml.dom.minidom.parse(enrichStats)
	
		runTag = xmlDom.getElementsByTagName('RunStats')[0]
		errTag = runTag.getElementsByTagName('ErrorRate')[0]
		r1Err = errTag.getElementsByTagName('DescriptiveStats')[0].getElementsByTagName('Mean')[0].childNodes[0].nodeValue
		r2Err = errTag.getElementsByTagName('DescriptiveStats')[1].getElementsByTagName('Mean')[0].childNodes[0].nodeValue

		runMetricArray['R1_MEAN_ERROR'] = round(float(r1Err), 2)
		runMetricArray['R2_MEAN_ERROR'] = round(float(r2Err), 2)
	else:
		fileMissing = 1

	hsMetricsList = os.listdir(align_folder)

	if len(hsMetricsList) > 0:

		recognizedPicardFormat = 0

		#Iterate through HSMetrics to withdraw GC Dropout values

		gcDropArray = []
		atDropArray = []

		for hsFile in hsMetricsList:
			if hsFile.endswith("HsMetrics.txt"):
				hsFileHandle = open(os.path.join(align_folder,hsFile), "r")
				for i, row in enumerate(hsFileHandle):
					if i == 6 and "BAIT_DESIGN_EFFICIENCY" in row:
						recognizedPicardFormat = 1
					if i == 7 and recognizedPicardFormat:
						gcDropArray.append(float(row.split("\t")[36]))
						atDropArray.append(float(row.split("\t")[35]))

		runMetricArray['GC_DROPOUT'] = round(numpy.mean(gcDropArray), 2)
		runMetricArray['AT_DROPOUT'] = round(numpy.mean(atDropArray), 2)
	else:
		fileMissing = 1

	if not fileMissing:
		conn = sqlite3.connect (metrics_db)
		conn.text_factory = str
		cur = conn.cursor()
		cur.execute("CREATE TABLE IF NOT EXISTS run_qc (ID VARCHAR, UNIFORMITY FLOAT, PCT_DUPLICATE FLOAT, TI_TV_RATIO FLOAT, TOTAL_PF_READS INTEGER, PCT_ALIGNED FLOAT, MEDIAN_FRAGMENT_LENGTH INTEGER, MEAN_COVERAGE_DEPTH FLOAT, PCT_Q30 FLOAT, PCT_ON_MANIFEST FLOAT, R1_MEAN_ERROR FLOAT, R2_MEAN_ERROR FLOAT, GC_DROPOUT FLOAT, AT_DROPOUT FLOAT, PRIMARY KEY (ID))")

		fieldNames = ['ID', 'UNIFORMITY', 'PCT_DUPLICATE', 'TI_TV_RATIO', 'TOTAL_PF_READS', 'PCT_ALIGNED', 'MEDIAN_FRAGMENT_LENGTH', 'MEAN_COVERAGE_DEPTH', 'PCT_Q30', 'PCT_ON_MANIFEST', 'R1_MEAN_ERROR', 'R2_MEAN_ERROR', 'GC_DROPOUT', 'AT_DROPOUT']

		fieldSQL = ",".join(fieldNames)
	
		valueList = []
	
		for field in fieldNames:
			valueList.append(str(runMetricArray[field]))

		valueSQL = ",".join(valueList)		
	
		sqlStatement = "INSERT OR REPLACE INTO run_qc (" + fieldSQL + ") VALUES (" + valueSQL + ")"
		cur.execute(sqlStatement)
		conn.commit()
		
def filterBAM(align_folder, sample_id, qscore, coverage_output_suffix, bam_suffix, samtoolsCMD):
	
	jobList = []
	
	originalBAM = glob.glob(os.path.join(align_folder,sample_id + "*" + bam_suffix + ".bam"))[0]
	processedBAM = originalBAM.replace(bam_suffix, coverage_output_suffix)
	processedLog = os.path.join(align_folder,sample_id + coverage_output_suffix + ".log")
			
	# If the processedBAM doesn't exist, use SAMTOOLS to dedup and filter
	if not os.path.isfile(os.path.join(align_folder,processedBAM)):
		print("\t\t" + sample_id + ".bam is being deduped and filtered by q" + str(qscore))
				
		# samtools view -bh -F d removes duplicates. while -qQSCORE filters MAPQ]		
		#job = subprocess.Popen(samtoolsCMD + " view -bh -F d -q" + str(qscore) + " -", stdout=open(processedBAM, 'w'), stderr=open(processedLog, 'w'), shell=True)
		command = " ".join([samtoolsCMD ,"view","-bh","-F","d","-q",str(qscore) ,"-o",processedBAM,originalBAM])
		#print(command)
		subprocess.call(command, stdout=open(processedBAM, 'w'), stderr=open(processedLog, 'w'), shell=True)
		
		# Append jobID and filehandlers to list
		#jobList.append(job)

	else:
		print("\t\t" + sample_id + ".bam has already been deduped and filtered.")
	
	# Clean-up once all the jobs are done
	#for job in jobList:
		#job.wait() 

def get_coverage(in_folder, analysis_folder, sample_id, qscore, coverage, genelist, bed, rscript, coverage_output_suffix, bedtoolsCMD, rscriptCMD):
	
	# Create output folder
	coverageOutFolder = os.path.join(analysis_folder,"coverage")
	geneCoverageOutFolder = os.path.join(coverageOutFolder,"gene_coverage")
	boxplotOutFolder = os.path.join(geneCoverageOutFolder,"boxplot")
	
	if not os.path.isdir(boxplotOutFolder):
		os.makedirs(boxplotOutFolder)
				
	bamFile = glob.glob(os.path.join(in_folder,sample_id + "*" + coverage_output_suffix + ".bam"))[0]
	histogramFile = os.path.join(coverageOutFolder,os.path.basename(bamFile).replace(".bam", ".hist"))
			
	#command = " ".join([bedtoolsCMD,"-abam",bamFile,"-b",bed,"-hist"]) # For 2.17
	command = " ".join([bedtoolsCMD,"coverage", "-a",bed,"-b",bamFile,"-hist"]) # For 2.26 (?) Website says -a and -b were flipped.
	#print command
	subprocess.call(command, stdout=open(histogramFile, 'w'), shell = True)
			
	[failed_exons_list, gene_exons_dict, exons_summary_csv] = calcMetrics(histogramFile, coverage, genelist)
			
	for geneCoverageFile in os.listdir(geneCoverageOutFolder):
		if geneCoverageFile.endswith(".txt"):
			geneCoveragePath = os.path.join(geneCoverageOutFolder, geneCoverageFile)
			geneCoverageLog = os.path.join(geneCoverageOutFolder, os.path.basename(geneCoverageFile).replace("txt","log"))

			# Use RScript to create coverage plots							
			command = " ".join([rscriptCMD, rscript, geneCoveragePath, boxplotOutFolder, str(coverage)])
			subprocess.call(command, stdout=open(geneCoverageLog, 'w'), stderr=subprocess.STDOUT, shell = True)

	return [failed_exons_list, gene_exons_dict, exons_summary_csv]
"""
=====================================
Helpers for get_coverage()
=====================================
def calcMetrics(inHistFile, depthThresh, geneListFile):
"""

def calcMetrics(inHistFile, depthThresh, geneListFile):
	
	failed_exons_list = []
	gene_exons_dict = {}
	
	interval = {}
	
	geneList = open(geneListFile)
	geneDict = {}
	
	for line in geneList:
		line = line.rstrip('\r\n')
		[geneName,accession] = line.split('\t')
		geneDict[accession] = geneName
		
	geneList.close()
	
	histogram = open(inHistFile, 'r')
	
	# Define output for 'all only' histogram
	allOnlyHist = open(inHistFile.replace('.hist','_all_only.hist'), 'w')
	
	# Read fields and coverage depth info for each region found in input histogram
	for line in histogram:		
		if 'all' in line:
			allOnlyHist.write(line)
		else:
			line = line.rstrip('\r\n')
			lineToken = line.split('\t')
			chr = lineToken[0]
			start = lineToken[1]
			end = lineToken[2]
			accession = '_'.join(lineToken[3].split('_')[0:2]) # Format NM_####_exon_#_chr#
			exon = '_'.join(lineToken[3].split('_')[2:4]) # Format NM_####_exon_#_chr#
			depth = int(lineToken[6])
			numAtDepth = int(lineToken[7])
			
			exonKey = '.'.join([chr,start,end,exon])
						
			if exonKey not in interval:
				interval[exonKey] = {}
				interval[exonKey]['depth'] = []
				interval[exonKey]['accession'] = accession
				interval[exonKey]['exon'] = exon
				interval[exonKey]['depthZero'] = 0
				interval[exonKey]['depthThresh'] = 0
				interval[exonKey]['sampleTotal'] = 0
				
			for i in range(numAtDepth):
				interval[exonKey]['depth'].append(depth)
			
			if depth == 0:
				interval[exonKey]['depthZero'] = interval[exonKey]['depthZero'] + numAtDepth
			

			if depth >= int(depthThresh):
				interval[exonKey]['depthThresh'] = interval[exonKey]['depthThresh']  + numAtDepth
			
			interval[exonKey]['sampleTotal'] = interval[exonKey]['sampleTotal'] + depth * numAtDepth
			
	histogram.close()
	
	# Define different output flies
	coverageTmp = open(inHistFile.replace('.hist','.tmp'), 'w')
	coverageTmp.write("chr\ttotal\tavg\tsample_total\tsample_avg\tmin\tq1\tmedian\tq3\tmax\tpctGT20\tpctEQ0\n")
	coverageSum = open(inHistFile.replace('.hist','.summary'), 'w')
	coverageSum.write("gene\trefSeq\texon\tchr\tstartPos\tendPos\ttotal\tavg\tsample.total\tsample.avg\tmin\tq1\tmedian\tq3\tmax\tpctGT20\tpctEQ0\n")

	reseqFailed = []
	reseqPassed = []
	
	geneSummary = {}
	
	for accession in geneDict:
		geneSummary[geneDict[accession]] = {}
	
	# Go through each exon region extracted from input histogram and write coverage stats to output
	for region in sorted(interval):
		numBase = len(interval[region]['depth'])
		total = interval[region]['sampleTotal']
		exon = interval[region]['exon']
		accession = interval[region]['accession']
		gene = geneDict[accession]
		
		q1idx = int(math.floor (0.25 * (numBase + 1)) - 1)
		q2idx = int(math.floor(0.5 * (numBase + 1)) - 1)
		q3idx = int(math.floor (0.75 * (numBase + 1)) - 1)

		min = str(interval[region]['depth'][0])
		q1 = str(interval[region]['depth'][q1idx])
		q2 = str(interval[region]['depth'][q2idx])
		q3 = str(interval[region]['depth'][q3idx])	
		max = str(interval[region]['depth'][-1])
		
		avg = int(math.floor(total  / numBase))
		#print '%.3f' % (float(int(interval[region]['depthThresh']) / int(numBase)))
		pctAboveThreshold = str(int(interval[region]['depthThresh']) * 100 / int(numBase))
		#pctAboveThreshold = str(numBase)
		pctZero = str(interval[region]['depthZero'] / numBase )
	
		[chr,start,end,exon] = region.split('.')

		if int(min) < int(depthThresh):
			filter = 'FAIL'
			reseqFailed.append(','.join([filter,gene,accession,exon,chr,start,end,min,pctAboveThreshold]))
			failed_exon = [filter,gene,accession,exon,chr,start,end,min,pctAboveThreshold,str(interval[region]['depthThresh']),str(numBase)]
			failed_exons_list.append(failed_exon)
			
			# Add to gene_exons_count 
			if gene not in gene_exons_dict:
				exon_summary = {}
				exon_summary['exons_total'] = 1
				exon_summary['exons_failed'] = 1
				exon_summary['exons_passed'] = 0
				exon_summary['exons_min_gt_20'] = 0
				exon_summary['exons_min_gt_50'] = 0
				exon_summary['exons_min_gt_100'] = 0
				
				gene_exons_dict[gene] = exon_summary
			else:
				gene_exons_dict[gene]['exons_total'] = gene_exons_dict[gene]['exons_total']  + 1
				gene_exons_dict[gene]['exons_failed'] = gene_exons_dict[gene]['exons_failed']  + 1
				
		else:
			filter ='PASS'
			reseqPassed.append(','.join([filter,gene,accession,exon,chr,start,end,min,pctAboveThreshold]))
			
			#print "PctAboveThreshold"
			#print pctAboveThreshold
			
			# Add to gene_exons_count 
			if gene not in gene_exons_dict:
				exon_summary = {}
				exon_summary['exons_total'] = 1
				exon_summary['exons_failed'] = 0
				exon_summary['exons_passed'] = 1
				exon_summary['exons_min_gt_20'] = 0
				exon_summary['exons_min_gt_50'] = 0
				exon_summary['exons_min_gt_100'] = 0
				
				gene_exons_dict[gene] = exon_summary
				
				if int(min) > 20:
					exon_summary['exons_min_gt_20'] = 1
				else:
					exon_summary['exons_min_gt_20'] = 0
					
				if int(min) > 50:
					exon_summary['exons_min_gt_50'] = 1
				else:
					exon_summary['exons_min_gt_50'] = 0
				
				if int(min) > 100:
					exon_summary['exons_min_gt_100'] = 1
				else:
					exon_summary['exons_min_gt_100'] = 0
			else:
				gene_exons_dict[gene]['exons_total'] = gene_exons_dict[gene]['exons_total']  + 1
				gene_exons_dict[gene]['exons_passed'] = gene_exons_dict[gene]['exons_passed']  + 1
				
				if int(min) > 20:
					exon_summary['exons_min_gt_20'] = exon_summary['exons_min_gt_20'] + 1
				if int(min) > 50:
					exon_summary['exons_min_gt_50'] = exon_summary['exons_min_gt_50'] + 1
				if int(min) > 100:
					exon_summary['exons_min_gt_100'] = exon_summary['exons_min_gt_100'] + 1
				

		coverageTmp.write(chr + ":" + start + "-" + end + "\t" + str(total) + "\t" + str(avg) + "\t" + str(total) + "\t" + str(avg) + "\t" + min + "\t" + q1 + "\t" + q2+ "\t" + q3 + "\t" + max + "\t" + pctAboveThreshold + "\t" + pctZero + "\n")
		coverageSum.write('\t'.join([gene ,accession,exon,chr,start,end,str(total),str(avg),str(total),str(avg),min,q1,q2,q3,max,pctAboveThreshold ,pctZero]) + "\n")
		
		# need to convert exon to exonNum - otherwise it will sort by alphabetical instead of numerical order
		#print exon
		exonNum = int(exon.split('_')[1])
		geneSummary[gene][exonNum] = '\t'.join([gene ,accession,exon,chr,start,end,str(total),str(avg),str(total),str(avg),min,q1,q2,q3,max,pctAboveThreshold ,pctZero])

	coverageTmp.close()
	coverageSum.close()
		
	# Write resequence summary - sorted by filter (FAIL first, PASS second)
	coverageReseq = open(inHistFile.replace('.hist','_resequence.csv'), 'w')
	exons_summary_csv = inHistFile.replace('.hist','_resequence.csv')
	coverageReseq.write("filter,gene,refSeq,exon,chr,startPos,endPos,min,pctGT20\n")
	
	for failed in reseqFailed:
		coverageReseq.write(failed + "\n")
	for passed in reseqPassed:
		coverageReseq.write(passed + "\n")

	coverageReseq.close()
	
	# Write gene-specific summary - sorted by exon instead of chr coordinates	
	
		# Open file handle for individual genes	
	coverageRootDir = os.path.dirname(inHistFile)
	geneFileHandleDict = {}
	
	for accession in geneDict:
		gene = geneDict[accession]
		
		if not os.path.isdir(os.path.join(coverageRootDir,'gene_coverage')):
			os.mkdir(os.path.join(coverageRootDir,'gene_coverage'))
		
		coverageBaseFile = os.path.basename(inHistFile)
		coverageGeneFile = coverageBaseFile.replace('.hist','_' + gene + '.txt')
		geneFileHandle = open(os.path.join(coverageRootDir, 'gene_coverage', coverageGeneFile), 'w')
		geneFileHandleDict[gene] = geneFileHandle
		geneFileHandleDict[gene].write("gene\trefSeq\texon\tchr\tstartPos\tendPos\ttotal\tavg\tsample.total\tsample.avg\tmin\tq1\tmedian\tq3\tmax\tpctGT20\tpctEQ0\n")

	for gene in geneSummary:
		for exon in sorted(geneSummary[gene]):
			geneFileHandleDict[gene].write( geneSummary[gene][exon] + "\n")

	for gene in geneFileHandleDict:
		geneFileHandleDict[gene].close()
		
	return [failed_exons_list, gene_exons_dict, exons_summary_csv]
		

def flatten_av(in_av, out_flat_av, av_info, av_geno):
	
	infoFieldsList = av_info.split(',')
	genoFieldsList = av_geno.split(',')
	
	#print in_av
	#print out_flat_av
	
	infile = open(in_av)
	outfile = open(out_flat_av, 'w')

	# Quickly jump to line 2 to read the file headers
	header = ""

	for i, line in enumerate(infile):
		
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')		

		if i == 0:
			header = line

			# Replace if these headers are kept
			#header = header.replace("ljb2_all", "SIFT\tPP2-HVar\tPP2-HVar-Category\tPP2-HDiv\tPP2-HDiv-Category\tLRT\tLRT-Category\tMT\tMT-Category\tMA\tMA-Category\tFATHMM\tGERP++\tPhyloP\tSiPhy")
			header = header.replace("Otherinfo", "Zygosity\tQUAL.AV\tDP.AV\tCHR\tCHRPOS\tdbSNPID\tREF\tOBS\tQUAL\tFILTER")
		if i == 1:
			infoCol = lineSplit[-3]
			genoHeader = lineSplit [-2]
			genoCol = lineSplit[-1]

			infoDict = parseAndFilterInfoValues(infoCol, infoFieldsList)
			genoDict = parseAndFilterGenoValues(genoHeader, genoCol, genoFieldsList)	
			
			for infoKey in sorted(infoDict.keys()):
				header = header + "\t" + infoKey

			for genoKey in sorted(genoDict.keys()):
				header = header + "\t" + genoKey

		if i > 1:
			break
	
	infile.close()
	
	outfile.write(header + "\n")	
	
	infile = open(in_av)
	
	# Print body
	for i, line in enumerate(infile):
		if i != 0:
			line = line.rstrip('\n\r')
			lineSplit = line.split('\t')
		
			# Join everything that doesn't need post-processing
			line = '\t'.join(lineSplit[0:-12])
			
			# Join the Misc info columns
			miscCols = '\t'.join(lineSplit[-12:-3])
			miscCols = miscCols + '\t'

			# Split INFO column into individual fields
			infoChunk = lineSplit[-3]
			infoCols = ""
			infoDict = parseAndFilterInfoValues(infoChunk, infoFieldsList)
			
			for infoKey in sorted(infoDict.keys()):
				infoCols = infoCols + infoDict[infoKey] + "\t"
			
			# Split GENO Columns into individual fields
			genoHeadChunk = lineSplit[-2]
			genoValChunk = lineSplit[-1]
			genoCols = ""
				
			genoDict = parseAndFilterGenoValues(genoHeadChunk, genoValChunk, genoFieldsList)

			for genoKey in sorted(genoDict.keys()):
				genoCols = genoCols + genoDict[genoKey] + "\t"
				
			genoCols = genoCols.replace("/",",") # This fixes the GT field (1/1) which becomes a date when uploaded to excel
			
			line = line + "\t" + miscCols + infoCols + genoCols + "\n"

			line = line.replace("\t\n","\n")
			
			outfile.write(line)

	infile.close()
	outfile.close()


"""
=====================================
Helpers for flatten_av()
=====================================
def parseAndFilterInfoValues(infoCol, interestedInfoFields):
def parseAndFilterGenoValues(genoHeaderCol, genoValCol, interestedGenoFields):
"""

def parseAndFilterInfoValues(infoCol, interestedInfoFields):
	
	infoDict = dict.fromkeys(interestedInfoFields, "NA") # Dictionary 
	
	infoArray= infoCol.split(';')

	for idx in enumerate(infoArray):
		num = idx[0]
		value = idx[1]
	
		if "=" in value:
			infoHeader = value.split('=')[0]
			infoVal = value.split('=')[1]
			
			if infoHeader in infoDict:
				infoDict[infoHeader] = infoVal
		
		"""
			"QD" in infoHeader: # quality score divided by total depth. PASS: QD > 2.0
			"BaseQRankSum" in infoHeader:
			"DP" in infoHeader: # Total depth of reads. PASS: DP > 20
			"HaplotypeScore"
			"ReadPosRankSum" in infoHeader: # Chance of it being FP cuz it's end of reads. ReadPosRankSum > -8
			"MQ" in infoHeader: # Mapping Quality. MQ > 30
			"VQSLOD" in infoHeader: # log odd ratios of it being a true variants. 
			"SB" in infoHeader: #Strand bias: higher = higher bias (Flag if value is >=0?)
			"FS" in infoHeader: #Fisher Strand: Indication of strand bias. FS < 60
		"""

	return infoDict

def parseAndFilterGenoValues(genoHeaderCol, genoValCol, interestedGenoFields):
	genoDict = dict.fromkeys(interestedGenoFields, "NA")
	genoHeader = genoHeaderCol.split(':')
	genoVal = genoValCol.split(':')

	if len(genoHeader) == len(genoVal):
		for idx in enumerate(genoHeader):
			num = idx[0]
			header = idx[1]
			if header in genoDict:
				genoDict[header] = genoVal[num]

			"""
			"GT" # 0/0  - Homo reference 0/1 - Hetero 1/1 Homo alternate
			"AD(Ref,Alt)" = genoArray[idx] # Depth per allele - REF depth / ALT depth
			"DP" # Total depth of reads - Flag if DP < 20
			"GQ" # Genotype quality. Phred-scale confidence that GT is true. Flag if < 30
			"SB" # Strand Bias. Flag as strand bias if value is >= 0
	
			"""
			
	# Added fix for empty DP
	alleleDepth = genoDict['AD'].split(',')
	alleleDepth = [int(i) for i in alleleDepth]
	genoDict['DP'] = str(sum(alleleDepth))
		
	return genoDict

def filter_av(av_in, av_out, genelist, bed, av_fields, variants_db, cnv_only_dict):
	
	#conn = sqlite3.connect('/home/miseq/miseq_pipe/resources/csv_db/variants.db')
	#print variants_db
	
	conn = sqlite3.connect(variants_db)
	conn.row_factory = sqlite3.Row
	
	avFieldsList = av_fields.split(',')
	bedDict = getBedDict(bed)
	geneList= getAccessionNumbers(genelist)

	infile = open(av_in)
	outfile = open(av_out, 'w')
	
	# Keeps track of which columns to list (based on avFieldsList)
	colIdxDict = {}
	colHeadingDict = {}
	
	# Read the first line of inFile and store its heading
	for i, line in enumerate(infile):
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		if i == 0:
			for j, heading in enumerate(lineSplit):
				colHeadingDict[heading] = j
				if heading in avFieldsList:
					colIdxDict[j] = heading
		if i > 0:
			break
	infile.close()
		
	# Print body
	infile = open(av_in)
	
	for i, line in enumerate(infile):
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		
		# Reset line to blank
		line = ""
		chr = lineSplit[colHeadingDict['Chr']].replace("chr","")
		chrCoord = lineSplit[colHeadingDict['Start']]
		ref = lineSplit[colHeadingDict['Ref']]
		alt = lineSplit[colHeadingDict['Alt']]
		geneName = lineSplit[colHeadingDict['Gene.refGene']].replace('(','').replace(')','').split(';')[0]
		av_transcripts = lineSplit[colHeadingDict['AAChange.refGene']].split(',')
		
		# Parse cDNA and AA from the correct transcript from AnnoVar
		
		av_cDNA = "c.?"
		av_AA = "p.?"
		
		for av_transcript in av_transcripts:
			av_transcript_fields = av_transcript.split(':')
			if len(av_transcript_fields) >= 3:
				av_accession = av_transcript_fields[1]
				if av_accession in geneList:
					av_cDNA = av_transcript_fields[3]
					if len(av_transcript_fields) >= 4:
						av_AA = av_transcript_fields[4]

		if i > 0:
			[accession, exon] = getTranscriptAndExonFromBedDict(chr,chrCoord,bedDict)
	
			if accession not in geneList:
				continue # skips the current line, go to the next iteration without processing it
		
			geneName = geneList[accession]
			
			if geneName in cnv_only_dict:
				continue

			if "c." not in av_cDNA or av_cDNA == "c.?":
				#print "Get HGCS From Mutalyzer"
				av_cDNA = get_hgvs_from_mutalyzer(chr, chrCoord, ref, alt, accession)

		# Print out heading / reformatted values from the flattened AnnoVar
		for j in sorted(colIdxDict.keys()):
				
			heading = colIdxDict[j]
			lineSplit[j] = lineSplit[j].replace(",", ";") #converting from TSB to CSV. Changing , to ;
			if heading == "Gene.refGene":
				if i == 0:
					line = line + "Gene.refGene.fixed" + ","
				else:
					line = line + geneName + ","
			elif heading == "AD":
				if i == 0:
					line = line + "AD.REF,AD.ALT" + ','
				else:
					adSplit = lineSplit[j].split(';')
					
					# Split the AD intp REF and ALT and then add it.
					line = line + adSplit[0] + ',' + adSplit[1] + ','
		
			# We skip Gene.refGene because it was replaced manually.
			# We skip AAChange.refGene because it will be manually processed in the next step.	
			#if heading == "AAChange.refGene" or heading == "Gene.refGene":
			if heading == "AAChange.refGene":
				continue
			else:
				line = line + lineSplit[j] + ','

		# Update Variant Fields Based on Variant DB and concatenate it to AnnoVar
		if i == 0:
			line = line + "AccessID.refGene,Exon.refGene,cDNAChange,AAChange,cDNAType,AAType,LatestCategory,ACMG,Omit,OmitComment,LatestDate,User"
		else:
		
			# Print variant key for debug
			#print("-".join([chr,chrCoord,ref,alt,accession]))

		
			#print (".".join([chr,chrCoord,ref,alt,accession]))

			variantsDict = getVariantsInfo_sqlite(conn,chr,chrCoord,ref,alt,accession)
			
			#print ",".join(chr,chrCoord,ref,alt)
			#print "Original"
			#print variantsDict

			# If no results found, try the reverse compliment
			if not variantsDict:
				#print "Original:" + ref + "," + alt
				ref, alt = complimentary_allele(ref, alt)				
				#print "revComp:" + ref + "," + alt
				variantsDict = getVariantsInfo_sqlite(conn,chr,chrCoord,ref,alt,accession)
				#print ",".join([chr,chrCoord,ref,alt])

			#print "revComp"
			#print variantsDict
				
			if variantsDict:
				line = concat_variantsDict_to_line(variantsDict, exon, accession, line)
				acmg = variantsDict['acmg']

			#else:
			#	if (exonDict[geneName]['strand'] == "+"):
			#		cdnaDict = chr2cds.find_cdna_fwd_strand_coords(exonDict, geneName, chrCoord)
			#	else: #if (exonDict[geneName]['strand'] == "-"):
			#		cdnaDict = chr2cds.find_cdna_rev_strand_coords(exonDict, geneName, chrCoord)
							
			#	strandDirection = exonDict[geneName]['strand']
				
				# Predict cDNA nomenclature
				### Replace this cDNA with the hgvs module
				### cDNA = predict_cdna(ref, alt, strandDirection, cdnaDict)
			else:
				acmg = 0
				line = line + ','.join([accession, exon, av_cDNA, av_AA, "-", "-", "Not Assessed", str(acmg),"N","-","-","-"])
					

			#reportCategory = acmgCat_to_reportCat(acmg)
			#line = line + ',' + reportCategory
			
		line = line + '\n'
			
		outfile.write(line)

	infile.close()
	outfile.close()
	
	conn.close()

"""
HGVS Helpers
def complimentary_allele(ref, alt):
def get_hgvs_from_mutalyzer(chromo, pos, ref, alt, accession):
def get_cnv_hgvs_from_mutalyzer(chromo, pos, ref, alt, accession):
"""

def complimentary_allele(ref, alt):

	allele_dict = {'G':'C', 'A':'T', 'T':'A', 'C':'G', '-':'-'}
	newRef = ""
	newAlt = ""

	for nuc in ref:
		newRef = newRef + allele_dict[nuc]

	for nuc in alt:
		newAlt = newAlt + allele_dict[nuc]

	return newRef, newAlt

def get_cnv_hgvs_from_mutalyzer(chromo, start_pos, end_pos, cnv_type, accession):
	
	variant = "c.?"
	
	#print "Attempt to get HGVS for " + ref + ">" + alt

	try: 
		URL =  'https://mutalyzer.nl/services/?wsdl'
		c = Client(URL, cache=None)
		o = c.service

		var_id = ""

		if "chr" not in str(chromo):
			chromo = "chr" + str(chromo)

		cnv_type_dict = {'deletion':'del', 'duplication':'dup', 'insertion':'ins'}

		var_id = chromo + ":g." + str(start_pos) + "_" + str(end_pos) + cnv_type_dict[cnv_type]

		r = o.numberConversion(build="hg19", variant=var_id)

		if r.string:
			for result in r.string:
				result_accession = result.split(":")[0].split(".")[0]
				result_variant = result.split(":")[1]
	
				if accession == result_accession:
					variant = result_variant
	except:
		pass
	
	return variant

def get_hgvs_from_mutalyzer(chromo, pos, ref, alt, accession):
	
	variant = "c.?"
	
	try: 
		URL =  'https://mutalyzer.nl/services/?wsdl'
		c = Client(URL, cache=None)
		o = c.service

		if "chr" not in str(chromo):
			chromo = "chr" + str(chromo)

		# DEL
		if len(ref) > 1:
			pos_end = int(pos) + (len(ref) - len(alt))
			var_id = chromo + ":g." + str(pos) + "_" + str(pos_end) + "del"
		# DEL of 1 nuc only
		elif alt == "-":
			var_id = chromo + ":g." + str(pos) + "del"
		# INS
		elif len(alt) > 1:
			pos_end = int(pos) + (len(alt) - len(ref))
			var_id = chromo + ":g." + str(pos) + "_" + str(pos_end) + "ins" + alt
		# INS of 1 nuc only
		elif ref == "-":
			pos_end = int(pos) + 1
			var_id = chromo + ":g." + str(pos) + "_" + str(pos_end) + "ins" + alt
		# SUB
		else:
			var_id = chromo + ":g." + str(pos) + ref + ">" + alt

		r = o.numberConversion(build="hg19", variant=var_id)


		if r.string:
			for result in r.string:
				result_accession = result.split(":")[0].split(".")[0]
				result_variant = result.split(":")[1]
	
				if accession == result_accession:
					variant = result_variant
	except:
		pass
	
	return variant

"""
=====================================
Helper modules for filter_av:
=====================================
def get_hgvs_from_mutalyzer(chromo, pos, ref, alt, accession)
def getBedDict(geneBedFile)
def getAccessionNumbers(geneListFile)
def concat_variantsDict_to_line(variantsDict, exon, accession, line)
def getVariantsInfo_sqlite(conn,chr,start_pos,ref,alt,transcript)
def getTranscriptAndExonFromBedDict(chr,chrPos,bedDict)
"""


def getBedDict(geneBedFile):
	bedDict= {}
	infile = open(geneBedFile)
	
	for line in infile:
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		chr = lineSplit[0].replace('chr','')
		chrStart = str(lineSplit[1])
		chrEnd = str(lineSplit[2])
		chrKey = chrStart + "-" + chrEnd
		exon = str(lineSplit[3].split('_')[3])
		accession = '_'.join(lineSplit[3].split('_')[0:2])
		
		if chr not in bedDict:
			bedDict[chr] = {}
		
		# For each chromosome, define a hash
		# The hash will contain chr start, chr end, and exon number
		bedDict[chr][chrKey] = {}
		bedDict[chr][chrKey]['exon'] = exon
		bedDict[chr][chrKey]['accession'] = accession
		
	infile.close()
	
	return bedDict
	
def getTranscriptAndExonFromBedDict(chr,chrPos,bedDict):
	exonNumber = "-"
	accession = "-"
	
	#for row in bedDict:
	#	print row
	
	# Need to check if chromosome and range is within dictionary. If not, then omit.
	if chr in bedDict:
		for chrRange in bedDict[chr]:
			chrStart = chrRange.split('-')[0]
			chrEnd = chrRange.split('-')[1]
		
			if chrPos > chrStart and chrPos <= chrEnd:
				exonNumber = bedDict[chr][chrRange]['exon']
				accession = bedDict[chr][chrRange]['accession']
			
	return [accession, exonNumber]
	
def getAccessionNumbers(geneListFile):
	geneList = {}
	infile = open(geneListFile)
	
	for line in infile:
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		geneName = lineSplit[0]
		accession = lineSplit[1]

		geneList[accession] = geneName

	infile.close()
	
	return geneList
	
def getVariantsInfo_sqlite(conn,chr,start_pos,ref,alt,transcript):
	
	acmgDict = {}
	
	acmgDict['Benign'] = 5
	acmgDict['Likely Benign'] = 4
	acmgDict['Predicted Benign'] = 4
	acmgDict['VUS'] = 3
	acmgDict['Uncertain Significance'] = 3
	acmgDict['Likely Pathogenic'] = 2
	acmgDict['Predicted Pathogenic'] = 2
	acmgDict['Pathogenic'] = 1

	variantsDict = {}
	cur = conn.cursor()
	
	sql = """SELECT * FROM variants WHERE chr='{0}' AND start_pos={1} AND ref='{2}' AND alt='{3}' AND transcript='{4}'""".format(chr,start_pos,ref,alt,transcript)
	cur.execute(sql)
	
	if cur:
		for row in cur:
			variantsDict = dict(itertools.izip(row.keys(), row))
	# This is newly added because ACMG is outdated in Access Database now
			if variantsDict['mutation_type'] in acmgDict:
				variantsDict['acmg'] = acmgDict[variantsDict['mutation_type']]
			else:
				variantsDict['acmg'] = 0
	
	return variantsDict
	
def concat_variantsDict_to_line(variantsDict, exon, accession, line):
	cDNA = variantsDict['cdna_change']
	aa_change = variantsDict['aa_change']
	cdna_type = variantsDict['cdna_type']
	aa_type = variantsDict['aa_type']
	mutation_type = variantsDict['mutation_type'].replace("Predicted Benign", "Pred. Benign").replace("Predicted Pathogenic", "Pred. Patho.")
	acmg = str(variantsDict['acmg'])
	omit = variantsDict['omit']
	omit_comment = variantsDict['omit_comment']
	#~ last_updated = variantsDict['date_last_updated'].strftime('%Y-%m-%d')
	last_updated = variantsDict['date_last_updated']
	user_updated  = variantsDict['user_updated']
				
	line = line + ','.join([accession,exon,cDNA,aa_change,cdna_type,aa_type,mutation_type,acmg,omit,omit_comment,last_updated,user_updated])
	
	return line
	
def report_qc_miseq(run_folder, run_id, sample_id, instrument):

	report_qc_dict = {}
	
	if instrument == "MiSeq":
		summary_file = run_folder + "/Data/Intensities/BaseCalls/Alignment/Summary.xml"
		summary_soup = BeautifulSoup(open(summary_file), "xml")
		summary = summary_soup.findAll("Summary")[0]	
		chipResultsSummary = summary.findAll("ChipResultsSummary")[0]
		report_qc_dict['yield_total']= "{0:.2f}".format(float(chipResultsSummary.findAll("yield")[0].string)/1000000000)
		report_qc_dict['run_date'] = run_id.split("_")[0]
		report_qc_dict['flowcell'] = run_id.split("_")[3].split("-")[1]
	else:
		report_qc_dict['yield_total'] = "NextSeq Yield"
		report_qc_dict['run_date'] = run_id.split("_")[0]
		report_qc_dict['flowcell'] = run_id.split("_")[3]
		
	my_qc = InteropDataset(run_folder)
	tile_qc = my_qc.TileMetrics()
	qual_qc = my_qc.QualityMetrics()
	q30List = qual_qc.read_qscore_results['q30']
	index_qc = my_qc.IndexMetrics()
	
	report_qc_dict['cluster_density'] = "{0:.2f}".format(float(tile_qc.mean_cluster_density_pf/1000))
	report_qc_dict['cluster_pf'] = "{0:.2f}".format(tile_qc.percent_pf_clusters)
	report_qc_dict['reads_pf']	= "{0:.2f}".format(index_qc.total_ix_reads_pf/1000000)
	report_qc_dict['pct_gt_q30'] = "{0:.2f}".format((sum(q30List) / len(q30List)))
	
	# Parse sample sheet
	
	report_qc_dict['investigator'] = "-"
	report_qc_dict['assay'] = "-"
	report_qc_dict['description'] = "-"
	report_qc_dict['experiment'] = "-"
	report_qc_dict['indexi5'] = "-"
	report_qc_dict['indexi7'] = "-"
	
	sampleSheetName = run_folder + "/SampleSheet.csv"
	sampleSheet = open(sampleSheetName, 'r')
	lines = sampleSheet.readlines()
	
	manifestDict = {}
	sampleSheetKeyDict = {}
	
	i5present = True #sometimes it's missing
	
	for line in enumerate(lines):
		# Check if i5 is present. Changes the parsing of run qc data
		if "Sample_ID" in line[1]:
			
			for idx,sampleKey in enumerate(line[1].rstrip('\r\n').split(',')):
				sampleSheetKeyDict[sampleKey] = idx
			
			if "I5_Index_ID" not in line[1]:
				i5present = False
				
		if "Investigator Name," in line[1]:
			report_qc_dict['investigator'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Assay Name," in line[1]:
			report_qc_dict['assay'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Experiment Name," in line[1]:
			report_qc_dict['experiment'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Description," in line[1]:
			report_qc_dict['description'] = line[1].split(',')[1].rstrip('\r\n')
		elif "[Manifests]" in line[1]:
			i = 1
			# Assertion: At least ONE manifest is present, else there's an error in the sample sheet
			manifestLine = lines[line[0]+i]
			
			while manifestLine not in ['\n', '\r\n']:
				manifestID = manifestLine.split(',')[0]
				manifestFile = manifestLine.split(',')[1].rstrip('\r\n')
				manifestDict[manifestID] = manifestFile
				i=i+1
				manifestLine = lines[line[0]+i]
				# This line error checks for empty lines filled with commas
				manifestLine = manifestLine.lstrip(',')
	
		elif sample_id == line[1].split(',')[0]:
			if (i5present):
				report_qc_dict['indexi5'] = line[1].split(',')[sampleSheetKeyDict['I5_Index_ID']] + ":" + line[1].split(',')[sampleSheetKeyDict['index2']]
				
			report_qc_dict['indexi7'] = line[1].split(',')[sampleSheetKeyDict['I7_Index_ID']] + ":" + line[1].split(',')[sampleSheetKeyDict['index']]
			
			if instrument == 'MiSeq':
				manifestID = line[1].split(',')[sampleSheetKeyDict['Manifest']]
				report_qc_dict['manifest'] = manifestDict[manifestID]
			else:
				report_qc_dict['manifest'] = line[1].rstrip('\r\n').split(',')[sampleSheetKeyDict['Description']]
	
	return report_qc_dict

def parse_cnv(sample_id, cnv_folder, cnv_panel_csv_file, bed_file, genelist_file):
	
	# Read the manifest bed into dictionary
	
	bed_dict = {}
	gene_dict = {}

	for line in open(genelist_file):
		line = line.rstrip('\r\n')
		[geneName,accession] = line.split('\t')
		gene_dict[geneName] = accession

	
	#print gene_dict

	for row in open(bed_file, 'r'):
		fields = row.split("\t")
		chr = fields[0]
		start = fields[1]
		end = fields[2]
		bed_key = start + "-" + end
		
		if chr not in bed_dict:
			bed_dict[chr] = {}
			
		if bed_key not in bed_dict[chr]:
			bed_dict[chr][bed_key] = ""
	
	cnv_dict = {}
	
	# Parse ExomeDepth results to cnv_dict
	exomedepth_filename = glob.glob(os.path.join(cnv_folder,"exomedepth",sample_id.replace("-",".")+"*.csv"))[0]
	exomedepth_file = csv.DictReader(open(exomedepth_filename, 'r'))

	for row in exomedepth_file:
		exomedepth_chr = "chr" + row['chromosome']
		
		exomedepth_in_manifest = 0
		
		# Check if call is within manifest
		if exomedepth_chr in bed_dict: 
			#print "ED Chr"
			#print exomedepth_chr
			#print bed_dict[exomedepth_chr]
			
			for bed_key in bed_dict[exomedepth_chr]:
			
				bed_start = int(bed_key.split("-")[0])
				bed_end = int(bed_key.split("-")[1])
				
				#print "BED_START " + str(bed_start) +  " ROW_START " + str(row["start"])
				#print "BED_END " + str(bed_end) + " ROW_END " + str(row["end"])		
				
				
				if int(row["start"]) >= bed_start and int(row["start"]) <= bed_end:
					exomedepth_in_manifest = 1
				elif int(row["end"]) >= bed_start and int(row["end"]) <= bed_end:
					exomedepth_in_manifest = 1
				elif int(row["start"]) <= bed_start and int(row["end"]) >= bed_end:
					exomedepth_in_manifest = 1

			if exomedepth_in_manifest:
				cnv_key = exomedepth_chr + ":" + row['start'] + "-" + row['end']
				cnv_dict[cnv_key] = {}
				cnv_dict[cnv_key]['chr'] = exomedepth_chr
				cnv_dict[cnv_key]['start'] = row['start']
				cnv_dict[cnv_key]['end'] = row['end']
				cnv_dict[cnv_key]['cnv_type'] = row['type']

				if None in row:
					if type(row[None]) is str:
						cnv_dict[cnv_key]['region'] = row['panels.hg19'] + ";" + row[None]
					else:
						cnv_dict[cnv_key]['region'] = row['panels.hg19'] + ";" + ";".join(row[None])
				else:
					cnv_dict[cnv_key]['region'] = row['panels.hg19']

				cnv_dict[cnv_key]['read_ratio'] = str(row['reads.ratio'])
				cnv_dict[cnv_key]['evidence'] = "Read Depth"

				# Get the first exon, and then split by _ to get the name of gene
				cnv_gene = cnv_dict[cnv_key]['region'].split(",")[0].split("_")[0]
				cnv_dict[cnv_key]['region'] = cnv_dict[cnv_key]['region'].replace(cnv_gene + "_", "Ex")
				cnv_accession = gene_dict[cnv_gene]
				cnv_hgvs = get_cnv_hgvs_from_mutalyzer(exomedepth_chr, row['start'], row['end'], row['type'], cnv_accession)
				cnv_dict[cnv_key]['hgvs'] = cnv_hgvs
				cnv_dict[cnv_key]['gene'] = cnv_gene
	manta_passed = 1
		
	try:		
		# Parse Manta results to cnv_dict
		manta_filename = glob.glob(os.path.join(cnv_folder,"manta",sample_id + "*manta.cnv.vcf"))[0]
		manta_file = open(manta_filename, 'r')
	except:
		manta_passed = 0

	if manta_passed:
		# Parse panel_exons.csv exon information
		#cnv_panel_csv = csv.DictReader(open(os.path.join(resource_folder, "bed_files", "cnv_panel_exons.csv"), 'r'))
		cnv_panel_csv = csv.DictReader(open(cnv_panel_csv_file, 'r'))

		exon_dict = {}
	
		for row in cnv_panel_csv:
			exon_chr = "chr" + row['chromosome']
			if exon_chr not in exon_dict:
				exon_dict[exon_chr] = {}
			
			exon_key = str(row['start']) + "-" + str(row['end'])
		
			if exon_key not in exon_dict[exon_chr]:
				exon_dict[exon_chr][exon_key] = row['name']

		# Goes through manta results
		for row in manta_file:
			if not row.startswith("#"):
				# PROCESS
				manta_chr = row.split("\t")[0]
				manta_start = int(row.split("\t")[1])
				manta_type = row.split("\t")[2].split(":")[0]
				manta_info = row.split("\t")[7]
				manta_end = int(manta_info.split(";")[0].replace("END=", ""))
			
				manta_region_list = []

				# Translate Manta_Type to CNV_Type
			
				if "MantaDEL" in manta_type:
					manta_type = "deletion"
				elif "MantaDUP" in manta_type:
					manta_type = "duplication"
				elif "MantaINS" in manta_type:
					manta_type = "insertion"
			
			
				# Check if call is within manifest
				manta_in_manifest = 0

				if manta_chr in bed_dict:
					for bed_key in bed_dict[manta_chr]:
						bed_start = int(bed_key.split("-")[0])
						bed_end = int(bed_key.split("-")[1])
					
						#print "Check BED"
						#print bed_start
						#print bed_end

						if int(manta_start) >= bed_start and int(manta_start) <= bed_end:
							manta_in_manifest = 1
						elif int(manta_end) >= bed_start and int(manta_end) <= bed_end:
							manta_in_manifest = 1
						elif int(manta_start) <= bed_start and int(manta_end) >= bed_end:
							manta_in_manifest = 1
				
				# Annotate manta results with panel_exons.csv exon information
				if manta_in_manifest:
					if manta_chr in exon_dict:
						for exon_key in exon_dict[chr]:

							#print "MANTA START:" + str(manta_start)
							#print "MANTA END:" + str(manta_end)

							#print "Check Regions"

							exon_start = int(exon_key.split("-")[0])
							exon_end = int(exon_key.split("-")[1])
							#print exon_start
							#print exon_end						
	
							if manta_start >= exon_start and manta_start <= exon_end:
								manta_region_list.append(exon_dict[chr][exon_key])
							elif manta_end >= exon_start and manta_end <= exon_end:
								manta_region_list.append(exon_dict[chr][exon_key])
							elif manta_start <= exon_start and manta_end >= exon_end:
								manta_region_list.append(exon_dict[chr][exon_key])

						#print manta_region_list

					# Compare manta region results with cnv region results
					for manta_region in manta_region_list:
						exomedepth_called = -1
	
						# Check if manta is in exomedepth cnv region
						for cnv_key in cnv_dict:
							if manta_region in cnv_dict[cnv_key]['region']:
								exomedepth_called = cnv_key
							
						if exomedepth_called >= 0:
							cnv_dict[exomedepth_called]['evidence'] = "Read Depth & Split Reads"
						else:
							cnv_key = manta_chr + ":" + str(manta_start) + "=" + str(manta_end)
							cnv_dict[cnv_key] = {}
							cnv_dict[cnv_key]['chr'] = manta_chr
							cnv_dict[cnv_key]['start'] = manta_start
							cnv_dict[cnv_key]['end'] = manta_end
							cnv_dict[cnv_key]['cnv_type'] = manta_type
							cnv_dict[cnv_key]['region'] = ",".join(manta_region)
							cnv_dict[cnv_key]['read_ratio'] = "N/A"
							cnv_dict[cnv_key]['evidence'] = "Split Reads"
	
	#else:
		#Do something to warn that manta failed. But proceed to generate reports anyway.

	return cnv_dict

def get_cnv_plots(sample_id, genelist, cnv_plot_folder):
	
	cnv_plots_list = []
	
	for line in open(genelist):
		line = line.rstrip('\r\n')
		[geneName,accession] = line.split('\t')
		
		if geneName:
			#print geneName
			if len(glob.glob(os.path.join(cnv_plot_folder, sample_id.replace("-", ".")  + "*_" + geneName + ".png"))) > 0:
				cnv_plot = glob.glob(os.path.join(cnv_plot_folder, sample_id.replace("-", ".")  + "*_" + geneName + ".png"))[0]
				#print ("cnv_plots_list.append " + cnv_plot)
				cnv_plots_list.append(cnv_plot)
	
	return cnv_plots_list
	
def generate_report(resource_folder, script_folder, report_folder, case_folder_root, flatten_av, cnv_dict, qc_dict, meta_dict, genes_exons_summary, failed_exons, exon_summary_csv, cnv_only_dict, seq_only_dict, alt_plot, base_plot, boxplots_list, cnv_plots_list, reanalysis):

	template_folder = os.path.join(script_folder,"html")
	report_template = ("report_template.html")
	env = Environment(loader=FileSystemLoader(template_folder))
	template = env.get_template(report_template)
	
	var_benign = []
	var_predicted_benign = []
	var_not_assessed = []
	var_confirm = []
	var_pop_benign = []
	
	variants_reported = 0
	variants_qc = 0
	cnvs_qc = len(cnv_dict)
	variants_incident = 0
	exons_qc = 0

	for gene in genes_exons_summary:
		#Add the CNV QC for Gene
		if gene not in seq_only_dict:
			exons_qc_gene = exons_qc + 1
		#Add the Exon QC for Gene
		exons_qc = exons_qc + int(genes_exons_summary[gene]["exons_total"])
	
		
	with open(flatten_av, 'rb') as f:
		in_av = csv.DictReader(f)
		for row in in_av:
	#		for gene in cnv_only_dict:
	#			if gene not in row:
			if row["Gene.refGene"] not in cnv_only_dict and row["LatestCategory"] <> '-':
				try:
					row["ALT_PCT"] = '{0:.2f}'.format(float(float(row["AD.ALT"])/(float(row["AD.REF"])+float(row["AD.ALT"]))))
				except:
					row["ALT_PCT"] = "-1"

				if row ['LatestCategory'] == 'Benign':
					var_benign.append(row)
					variants_incident = variants_incident + 1
				elif row['LatestCategory'] == 'Predicted Benign' or row['LatestCategory'] == 'Likely Benign':
					var_predicted_benign.append(row)
					variants_qc = variants_qc + 1
				elif row['LatestCategory'] == 'Not Assessed':
					if row['ExAC_nontcga_ALL'] != '.' and  float(row['ExAC_nontcga_ALL']) > 0.05:
						row["LatestCategory"] = 'Benign (ExAC)'
						var_pop_benign.append(row)
						variants_incident = variants_incident + 1
					else:
						var_not_assessed.append(row)
						variants_qc = variants_qc + 1
				else:
					var_confirm.append(row)
					variants_qc = variants_qc + 1
							
				variants_reported = variants_reported + 1
				
	meta_dict['variants_reported'] = variants_reported
	meta_dict['variants_qc'] = variants_qc
	meta_dict['cnvs_qc'] = cnvs_qc
	meta_dict['variants_incident'] = variants_incident
	meta_dict['exons_qc'] = exons_qc # exons in gene summary
								
	template_vars = {"meta_dict": meta_dict,
				"variant_confirm": var_confirm,
				"variant_not_assessed": var_not_assessed,
				"variant_predicted_benign": var_predicted_benign,
				"variant_pop_benign": var_pop_benign,
				"variant_benign": var_benign,
				"genes_exons_summary": genes_exons_summary,
				"failed_exons": failed_exons,
				"cnv_dict":cnv_dict,
				"qc_dict": qc_dict,
				"alt_plot": alt_plot,
				"base_plot": base_plot,
				"boxplots_list": boxplots_list,
				"cnv_plots_list": cnv_plots_list,
			}
					
	html_out = template.render(template_vars)

	meta_dict['sample_id']
	report_css = os.path.join(script_folder,"html","report_css.css")
	
	report_pdf = os.path.join(report_folder, meta_dict['sample_id'] + "_NGS_Report_L1.pdf")
	if reanalysis == 1:
		report_pdf = os.path.join(report_folder, meta_dict['sample_id'] + "_Reanalysis_Panel_NGS_Report_L1.pdf")

	HTML(string=html_out).write_pdf(report_pdf, stylesheets=[report_css])
	
	# Skip allocating copy of case report if root is not defined at the beginning
	#print case_folder_root
	if case_folder_root: 
		# Add in report code - need to parse it from analysis_type.csv - skip if report is empty from sheet or dictionary.
		case_folder = ""
		
		if reanalysis == 1:
			case_folder = os.path.join(case_folder_root, "REANALYSIS_REPORTS", meta_dict['sample_id'])
		else:
			#Case Folder Organization
			if meta_dict['sample_id'].startswith("B16") or meta_dict['sample_id'].startswith("B17"):
				case_folder = os.path.join(case_folder_root, "BREAST", meta_dict['sample_id'])
			elif meta_dict['sample_id'].startswith("C16") or meta_dict['sample_id'].startswith("C17"):
				case_folder = os.path.join(case_folder_root, "COLON", meta_dict['sample_id'])
			elif meta_dict['sample_id'].startswith("RS"):
				case_folder = os.path.join(case_folder_root, "RESEARCH", meta_dict['sample_id'])
			elif meta_dict['sample_id'].startswith("PT"):
				case_folder = os.path.join(case_folder_root, "PROFICIENCY TESTING", meta_dict['sample_id'])
			else:
				case_folder = os.path.join(case_folder_root, "UNSORTED_NGS_REPORTS", meta_dict['sample_id'])
		
		case_intermediates = os.path.join(case_folder, "NGS_Analysis")

		if not os.path.exists(case_intermediates):
			os.makedirs(case_intermediates)
		try:
			# Report PDF
			#---------------
			shutil.copy(report_pdf, case_folder)
			# Variant File
			#---------------
			shutil.copy(flatten_av, case_intermediates)
			#  QC File
			#----------
			qc_file = open(os.path.join(case_intermediates, meta_dict['sample_id'] + "_qc_stats.csv"), 'w')
			qc_file.write (",".join(qc_dict.keys()) + "\n")
			qc_value = []
			for key in qc_dict.keys():
				qc_value.append(qc_dict[key])
			qc_file.write(",".join(qc_value))
			qc_file.close()
			# Exon summary file
			#-----------------------
			shutil.copy(exon_summary_csv, case_intermediates)
			# CNV summary file
			#----------------------
			if len(cnv_dict) >= 1:
				cnv_file = open(os.path.join(case_intermediates, meta_dict['sample_id'] + "_cnv.csv"), 'w')

				cnv_header = ['gene','region','chr','start','end','cnv_type','hgvs','read_ratio','evidence']
				cnv_file.write(",".join(cnv_header) + "\n")

				for cnv_key in cnv_dict:
					cnv_value = []
					for key in cnv_header:
						cnv_value.append(cnv_dict[cnv_key][key])			
					cnv_file.write(",".join(cnv_value) + "\n")
				cnv_file.close()
			#CNV gene QC file
			#--------------------
			cnv_qc = open(os.path.join(case_intermediates, meta_dict['sample_id'] + "_cnvqc.csv"), 'w')
			for gene in genes_exons_summary:
				if gene not in seq_only_dict:
					if genes_exons_summary[gene]["exons_failed"] < 3:
						cnv_qc_gene = gene + "_CNV,Passed,Completed" #Fail if the gene has 3 or more exons that are failed.
					else:
						cnv_qc_gene = gene + "_CNV,Failed,QC Blood"
					cnv_qc.write(cnv_qc_gene + "\n")
			cnv_qc.close()
			#Helper scripts
			#-----------------
			shutil.copy(os.path.join(resource_folder,"tools","custom","ngs_toolkit","helpers","upload_ngs_results_to_access.vbs"), case_intermediates)
			shutil.copy(os.path.join(resource_folder,"tools","custom","ngs_toolkit","helpers","Click_me_to_upload_NGS_results_to_database.lnk"), case_folder)

		except:
			pass
main()
