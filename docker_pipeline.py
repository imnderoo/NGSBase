#!/usr/bin/python

import os, sys, subprocess
import shutil
import argparse #AW: Package for parsing arguments
import datetime, time
import re, glob

def main():
	#AW: This block of code is easy to reuse for any script that take in arguments
	#----
	parser = argparse.ArgumentParser(description = 'Script to read VCF into DB for ALBI project')
	parser.add_argument('in_folder', help = 'Full path to folder containing input FASTQs [or input file]')
	parser.add_argument('roi_bed_path', help = 'Full path to ROI BED - Should be in resources')
	parser.add_argument('out_folder_path', help = 'Full path to output folder')
	parser.add_argument('--prefix', help = 'Prefix to add to front of sample name. Ex: "RUNID_Sample.*"') 
	parser.add_argument('--start', help = 'Stage to start in [fastq|bam] DEFAULT: fastq') # AW: Optional argument: --start option
	parser.add_argument('--end', help = 'Stage to start in [bam|vcf] DEFAULT: vcf') 
	parser.add_argument('--intermediates', help = 'Keep Intermediates [Y|N] DEFAULT: N') 
	#parser.add_argument('--intermediates', help = 'Keep Intermediates [Y|N] DEFAULT: N') 

	# TO-DO: Proper handling for starting in fastq/bam or ending in bam/vcf

	#parser.add_argument('--in_file1', help = 'Name of input file 1') # AW: Each line defines the input parameter and the name to use


	#Ex: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder 
	#Ex2: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder --end bam
	#Ex3: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder 
	#Ex4: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder --start bam --prefix 160428_M01234_FLOWCE

	args = parser.parse_args() # AW: This line is needed for arguments to be recognized
	#----

	if args.start:
		print("Pipeline starting at " + args.start)
	else:
		print("Pipeline starting at fastq")

	if args.end:
		print("Pipeline ending at " + args.end)
	else:	
		print("Pipeline ending at vcf")
		
	if args.prefix:
		print("Pipeline prefix: " + args.prefix)

	if not os.path.exists(args.out_folder_path):
		os.makedirs(args.out_folder_path)

	AbsPath_to_BedFile = args.roi_bed_path # AW: Note the usage of args.parameter_namE

	# Look in the input folder and extract R1 FASTQ as set of input file.
	input_file_set = []

	for fastq in os.listdir(args.in_folder):
		#print fastq
		if "_R1_" in fastq and not "Undetermined" in fastq:
			input_file_set.append(fastq)
			
	#print input_file_set

	# Overwite default if input start is set to bam
	if args.start == "bam":
		input_file_set = []
		for bam in os.listdir(args.in_folder):
			if ".bam" in bam:
				input_file_set.append(bam)

	# SET GLOBAL VARIABLES (ex. PATHS TO RESOURCES)
	AbsPathToResources=os.path.abspath('/media/sf_resources')
	GenomeFile = "ucsc.hg19.fa" # AW genome

	AbsPathToAnnovDBFld=AbsPathToResources + '/annovar_DBs/'
	AbsPathToPubDBs=AbsPathToResources+'/gatk_DBs/'
	AbsPathToGenome=AbsPathToResources + '/genome/' + GenomeFile

	temp_dir = "/media/sf_temp"
	
	#GET ABSOLUTE PATHS OF INPUT 
	args.out_folder_path = os.path.abspath(args.out_folder_path)
	args.in_folder = os.path.abspath(args.in_folder)
	args.roi_bed_path = os.path.abspath(args.roi_bed_path)
	OutputFolder = args.out_folder_path + "/"
	

	#DOCKER VOLUME MOUNTS

	if not os.path.isdir(args.out_folder_path):
		os.makedirs(args.out_folder_path)
	
	#dockerVolume='-v /media/:/media/'
	dockerVolume='-v ' + AbsPathToResources + ':' + AbsPathToResources
	dockerVolume=dockerVolume+' '+'-v '+args.out_folder_path+':'+args.out_folder_path
	dockerVolume=dockerVolume+' '+'-v '+args.in_folder+':'+args.in_folder
	dockerVolume=dockerVolume+' '+'-v '+temp_dir+':'+temp_dir # This is needed to copy output to a local scratch for MANTA CNV
	
	#dockerVolume=dockerVolume+' '+'-v '+os.path.dirname(args.roi_bed_path)+':'+os.path.dirname(args.roi_bed_path)

	#DOCKER COMMANDS
	bwaCMD='docker run --rm '+dockerVolume+' nderoo324/bwa'
	samtoolsCMD='docker run --rm '+dockerVolume+' nderoo324/samtools'
	picardCMD='docker run --rm '+dockerVolume+' nderoo324/picard'
	annovarCMD='docker run --rm '+dockerVolume+' nderoo324/av'
	gatkCMD='docker run --rm '+dockerVolume+' nderoo324/gatk GenomeAnalysisTK'
	exomeDepthCMD = 'docker run --rm '+dockerVolume+' nderoo324/exomedepth'
	mantaCMD = 'docker run --rm '+dockerVolume+' nderoo324/manta'

	# DEBUG OPTIONS - 1 == TRUE, 0 == FALSE	
	run_bwa = 1
	run_samtools_merge = 1
	run_picard = 1
	run_gatk = 1
	run_annovar = 1
	run_cnv = 0
	run_manta = 0

	# START OF VCF CALLING SCRIPT
	start_time = time.asctime( time.localtime(time.time()) )

	sample_dict = {}
	elapsed_time_dict = {}

	for inFile in input_file_set:
		######## The following quotation marks should be filled with the FastQ names
		#FastQ_R1_FullPath = args.fastq_R1_path
		#FastQ_R2_FullPath = args.fastq_R2_path
		FastQ_R1_FullPath = args.in_folder + "/" + inFile
		FastQ_R2_FullPath = args.in_folder + "/" + inFile.replace("R1", "R2")
		
		FastQ_R1 = inFile #AW: os.path.basename leaves us with just the file name
		FastQ_R2 = inFile.replace("R1", "R2")

		####### SAMPLE specific info ############
		#Sample_Number=file1[0:8]
		#Sample_Number = "_".join(FastQ_R1.split("_")[:2]) #If we want the S# straign from illumina MiSeq. It splits by _, takes first 2, joins by _
		#Sample_Number=FastQ_R1.split("_")[0] #If we want the S# straign from illumina MiSeq
		#print ("\nSample Number: ",Sample_Number)

		### Parse FASTQ file using REGEX ###

		#print inFile

		matches = re.match( r'(.*?)_L(.*?)_R.*', inFile)
		Sample_Number = matches.group(1)
		lane_id = "L" + matches.group(2)
		
		Platform="Illumina"
		Library="TruSightOne"
		PU_name="1"
		
		# Add prefix if defined
		if args.prefix:
			SampleLaneIDsegment = args.prefix + "_" + Sample_Number + "-" + lane_id
			SampleIDsegment = args.prefix + "_" + Sample_Number
		else:
			SampleLaneIDsegment = Sample_Number + "-" + lane_id
			SampleIDsegment = Sample_Number
			
		sample_start_time = datetime.datetime.now()
			
		#print "Output: " + OutputFolder

		OutputBamSrtFile = OutputFolder + SampleLaneIDsegment + "-aln-pe-sorted.bam"
		OutputBamSrtDeDupFile = OutputFolder+SampleIDsegment+"-aln-pe-sorted-merged-dedup.bam"
		
		if os.path.exists(OutputBamSrtDeDupFile):
			print ("******* SKIPPING BWA, MERGE and SORT on:", SampleIDsegment, ". SORTED-DEDUP BAM ALREADY EXISTS *******")
			run_bwa = 0	
			run_samtools_merge = 0
			run_picard = 0

		#Define and vopen log file for rw
		log_date = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
		OutputLog = OutputFolder+SampleIDsegment+"-pipeline-"+log_date+".log"

		if args.start:
			log_write(OutputLog, "Pipeline starting at " + args.start + "\n")
		else:
			log_write(OutputLog, "Pipeline starting at fastq\n") 
		if args.end:
			log_write(OutputLog, "Pipeline ending at " + args.end + "\n")
		else:	
			log_write(OutputLog, "Pipeline ending at vcf\n")	
		if args.prefix:
			log_write(OutputLog, "Pipeline prefix: " + args.prefix + "\n")

		#	print matches.group()
		#	print (sample_id, lane_id)
			


		## BWA ############
		if run_bwa:
			if args.start != 'bam':
				print ("********************* RUNNING BWA on:",FastQ_R1," and ",FastQ_R2," *******************")
				log_write (OutputLog, '********************* RUNNING BWA on:",FastQ_R1," and ",FastQ_R2," *******************\n')
				# For newest samtools - pipe is tough on docker / space			
				#Intermediate files
				tmpSAM = OutputFolder + SampleLaneIDsegment + ".tmp.sam"
				tmpBAM = OutputFolder + SampleLaneIDsegment + ".tmp.bam"
				#bwa and samtools
				subprocess.call([bwaCMD+" mem -t 8 -a \
				-R '@RG\\tID:"+SampleLaneIDsegment+"\\tSM:"+Sample_Number+"\\tPL:"+Platform+"\\tLB:"+Library+"\\tPU:"+PU_name+"' \
				"+AbsPathToGenome+" "+FastQ_R1_FullPath+" "+FastQ_R2_FullPath+" > "+tmpSAM], shell=True)
				subprocess.call([samtoolsCMD+" view -Shu -@ 5 "+tmpSAM+" > "+tmpBAM], shell=True) 
				subprocess.call([samtoolsCMD+" sort -@ 5 -o "+OutputBamSrtFile+" "+tmpBAM], shell=True)		
				# Clean up tmp files
				os.remove(tmpSAM)
				os.remove(tmpBAM)
				print ("********************* BWA FINISHED. OUTPUT:",OutputBamSrtFile," *******************")

		if SampleIDsegment not in sample_dict:
			sample_dict[SampleIDsegment] = []
			elapsed_time_dict[SampleIDsegment] = datetime.datetime.now() - sample_start_time
		else:
			elapsed_time_dict[SampleIDsegment] = elapsed_time_dict[SampleIDsegment] + (datetime.datetime.now() - sample_start_time)		
			
		sample_dict[SampleIDsegment].append(SampleLaneIDsegment)		
		#print sample_dict
		
	## Per Sample Basis ##
	for SampleIDsegment in sample_dict:
		
		sample_start_time = datetime.datetime.now()
		
		OutputBamMrgFile = OutputFolder + SampleIDsegment + "-aln-pe-sorted-merged.bam"
		OutputBamSrtDeDupFile = OutputFolder+SampleIDsegment+"-aln-pe-sorted-merged-dedup.bam"
		OutputBamSrtDeDupRecalFile = OutputFolder+SampleIDsegment+"-aln-pe-sorted-merged-dedup-recal.bam"
		RawVariantsVCF = OutputFolder+SampleIDsegment+"_raw_variants.vcf"

		#Define and vopen log file for rw
		log_date = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
		OutputLog = OutputFolder+SampleIDsegment+"-pipeline-"+log_date+".log"
		
		#print sample_dict[sample]
		# Merge the multiple lane files (if they exist)
		
		#if os.path.exists(OutputBamSrtDeDupFile):
		#	print ("****** SKIPPING MERGE AND DEDUP on "+SampleIDsegment+". -merged-dedup.bam already exists. *******")
		#	run_samtools_merge = 0
		#	run_picard = 0

		if run_samtools_merge:
			if len(sample_dict[SampleIDsegment]) == 1:
				SampleLaneIDsegment = sample_dict[SampleIDsegment][0]
				PartBamSrtFile = OutputFolder + SampleLaneIDsegment + "-aln-pe-sorted.bam"
				
				#print ("os.rename(" + PartBamSrtFile + ", " + OutputBamMrgFile + ")")
				
				if os.path.isfile(PartBamSrtFile):
					os.rename(PartBamSrtFile, OutputBamMrgFile)
				else:
					print("Error: " + PartBamSrtFile + " does not exist.")
			else:
				partBAM = []
				
				for part in sample_dict[SampleIDsegment]:
					SampleLaneIDsegment = part
					PartBamSrtFile = OutputFolder + SampleLaneIDsegment + "-aln-pe-sorted.bam"
					partBAM.append(PartBamSrtFile)
				
				subprocess.call([samtoolsCMD+" merge "+OutputBamMrgFile+" "+" ".join(partBAM)],shell=True)
			
		## Picard MarkDuplicates
		if run_picard:
			print ("******************* RUNNING Picard MarkDuplicate on "+SampleIDsegment+" *****************")
			log_write (OutputLog, "******************* RUNNING Picard MarkDuplicate on "+SampleIDsegment+" *****************\n")

			subprocess.call(\
			[picardCMD+" MarkDuplicates \
			INPUT="+OutputBamMrgFile+"  \
			OUTPUT="+OutputBamSrtDeDupFile+"  \
			METRICS_FILE="+OutputFolder+"/"+SampleIDsegment+"_PCR_duplicates  \
			CREATE_INDEX=true \
			REMOVE_DUPLICATES=true 2>&1 | tee -a "+OutputLog]
			,shell=True)
		
		if "bam" not in str(args.end):
			
			## IndelRealigner - No longer necessary with Haplotype 
			if run_gatk:
				############### Base Recalibration ###########################
				print ("***************** RUNNING GATK BaseRecalibrator on "+SampleIDsegment+" *********************")
				log_write (OutputLog, "***************** RUNNING GATK BaseRecalibrator on "+SampleIDsegment+" *********************\n")
				subprocess.call(
				gatkCMD+" -T BaseRecalibrator \
				-R "+AbsPathToGenome+" \
				-I "+OutputBamSrtDeDupFile+"  \
				-L "+args.roi_bed_path+" \
				-knownSites "+AbsPathToPubDBs+"1000Genomes/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
				-knownSites "+AbsPathToPubDBs+"1000Genomes/1000G_phase1.indels.hg19.sites.vcf \
				-knownSites "+AbsPathToPubDBs+"dbSNP/dbsnp_138.hg19.vcf \
				-o  "+OutputFolder+"/"+SampleIDsegment+"_recal_data.grp"+" 2>&1 | tee -a "+OutputLog\
				,shell=True) 
				
				# PrintReads -Step 2 of Base Recalibration # Added --filter_bases_not_stored
				print ("**************** RUNNING GATK PrintReads on "+SampleIDsegment+" *****************")
				log_write (OutputLog, "**************** RUNNING GATK PrintReads on "+SampleIDsegment+" *****************\n")
				subprocess.call( \
				gatkCMD+" -T PrintReads \
				-R "+AbsPathToGenome+" \
				-I "+OutputBamSrtDeDupFile+"  \
				-L "+args.roi_bed_path+ " \
				--filter_bases_not_stored \
				-BQSR "+OutputFolder+"/"+SampleIDsegment+"_recal_data.grp  \
				-o  "+OutputBamSrtDeDupRecalFile+" 2>&1 | tee -a "+OutputLog \
				,shell=True)
				
				#################  HaplotypeCaller ###########################
				print ("**************** RUNNING GATK HaplotypeCaller on "+SampleIDsegment+" ********************")
				log_write (OutputLog, "**************** RUNNING GATK HaplotypeCaller on "+SampleIDsegment+" ********************\n")
				subprocess.call( \
				gatkCMD+ " -T HaplotypeCaller \
				-R "+AbsPathToGenome+" \
				-I "+OutputBamSrtDeDupRecalFile+" \
				-L "+args.roi_bed_path+" \
				--genotyping_mode DISCOVERY \
				-stand_emit_conf 10 \
				-stand_call_conf 30 \
				-mbq 20 \
				-o "+RawVariantsVCF+" 2>&1 | tee -a "+OutputLog \
				,shell=True)
			
			if run_annovar:
				################ ANNOVAR #################################
				print ("**************** RUNNING ANNOVAR on "+SampleIDsegment+" ********************\n")
				log_write (OutputLog, "**************** RUNNING ANNOVAR on "+SampleIDsegment+" ********************\n")
				#VCF file as input to ANNOVAR table_annovar.pl
				# -protocol refGene,avsnp144,exac03nontcga,1000g2015aug_all,esp6500siv2_all,clinvar_20160302,nci60,cosmic70,ljb26_all \
				
				# Comphrehensive - For research or later purposes
				"""
				subprocess.call(\
				annovarCMD+" table_annovar.pl "+RawVariantsVCF+" "+AbsPathToAnnovDBFld+" -buildver hg19 \
				-protocol refGene,avsnp144,exac03nontcga,1000g2015aug_all,esp6500siv2_all,clinvar_20160302,ljb26_all \
				-operation g,f,f,f,f,f,f \
				-argument -hgvs,,,,,, \
				-out "+OutputFolder+"/"+SampleIDsegment+"_Annov_out  \
				-remove \
				-otherinfo \
				-nastring . \
				-vcfinput "+" 2>&1 | tee -a "+OutputLog \
				,shell=True)		
				"""

				subprocess.call(\
				annovarCMD+" table_annovar.pl "+RawVariantsVCF+" "+AbsPathToAnnovDBFld+" -buildver hg19 \
				-protocol refGene,exac03nontcga \
				-operation g,f \
				-argument -hgvs, \
				-out "+OutputFolder+"/"+SampleIDsegment+"_Annov_out  \
				-remove \
				-otherinfo \
				-nastring . \
				-vcfinput "+" 2>&1 | tee -a "+OutputLog \
				,shell=True)
			
		elapsed_time_dict[SampleIDsegment] = elapsed_time_dict[SampleIDsegment] + (datetime.datetime.now() - sample_start_time)
		
		log_write (OutputLog, "Time Elapsed: " + str(elapsed_time_dict[SampleIDsegment]) + "\n")

	### CNV CALLING SCRIPTS (ExomeDepth and Manta) ###
	
	if run_cnv:
		exomedepth_bam_list = glob.glob(OutputFolder+"/*-aln-pe-sorted-merged-dedup.bam") #recal is the trimmed down version of BAM - may miss out on some intergenic reads
		# recal might be needed for MANTA (to reduce the file count, but dedup is best for exomeDepth to keep everything
		
		# Manta requires local copy of the input files to avoid OS/Filesystem delays!!!
		manta_bam_list = glob.glob(OutputFolder+"/*-aln-pe-sorted-merged-dedup-recal.bam") #recal is the trimmed down version of BAM - may miss out on some intergenic reads
		

		for original_bam_file in manta_bam_list:
			temp_bam = original_bam_file.replace(OutputFolder.rstrip("/"), temp_dir) 
			temp_bai = original_bam_file.replace(OutputFolder.rstrip("/"), temp_dir).replace("bam", "bai")
			
			if run_manta:
				if not os.path.isfile(temp_bam):
					shutil.copy(original_bam_file, temp_dir)
				
				if not os.path.isfile(temp_bai):
					shutil.copy(original_bam_file.replace("bam","bai"), temp_dir)	
			
		manta_bam_list = glob.glob(temp_dir+"/*-aln-pe-sorted-merged-dedup-recal.bam")
		
		bam_suffix = "-aln-"
		cnv_panel_csv = os.path.join(AbsPathToResources, "bed_files", "cnv_panel_exons.csv")
		call_cnv(exomedepth_bam_list, manta_bam_list, OutputFolder, cnv_panel_csv, AbsPathToGenome, exomeDepthCMD, mantaCMD, bam_suffix, temp_dir, run_manta)
		
	# CLEAN OUT INTERMEDIATES
	if not args.intermediates:
		try:
				for f in os.listdir (OutputFolder):
					if  re.search("-aln-pe-sorted.ba", f) or \
					re.search("-aln-pe-sorted-merged.ba", f) or \
					re.search("PCR_duplicates", f) or \
					re.search("recal_data.grp", f) :
						os.remove(os.path.join(OutputFolder, f))
		except:
			pass
		
	end_time = time.asctime( time.localtime(time.time()) )
	
	# END OF SCRIPT
	print ("\nStart time: " + start_time + "\n")
	print ("End time: " + end_time + "\n")
	
def call_cnv(exomedepth_bam_list, manta_bam_list, out_dir, panel_csv, genome, exomeDepthCMD, mantaCMD, bam_suffix, temp_dir, run_manta):

	cnv_out = os.path.join(out_dir,"cnv")
	
	if not os.path.isdir(cnv_out):
		os.makedirs(cnv_out)
		
	# Create bamlist file
	bam_list_file = os.path.join(cnv_out, "exomedepth-bamlist.txt")

	with (open(bam_list_file, 'w')) as f:
		for in_bam in exomedepth_bam_list:
			f.write(in_bam + "\n")
	
	# ExomeDepth - This will run all bams at once.
	print ("**************** RUNNING ExomeDepth to call CNV ********************\n")
	#log_write (OutputLog, "**************** RUNNING ExomeDepth to call CNV ********************\n")
	
	exomedepth_out = os.path.join(cnv_out,"exomedepth")
	
	if not os.path.isdir(exomedepth_out):
		os.makedirs(exomedepth_out)
	
	subprocess.call(" ".join([exomeDepthCMD, \
		"-b", bam_list_file, \
		"-o", exomedepth_out, \
		"-p", panel_csv, \
		"-s", "0.3"]), shell=True)
	
	# Manta - Call CNV for each bam in the bam list.
	
	# Create bamlist file
	bam_list_file = os.path.join(cnv_out, "manta-bamlist.txt")

	with (open(bam_list_file, 'w')) as f:
		for in_bam in manta_bam_list:
			f.write(in_bam + "\n")
			
	if run_manta:
		print ("**************** RUNNING Manta to call CNV ********************\n")
		#log_write (OutputLog, "**************** RUNNING Manta to call CNV ********************\n")
	
		#manta_out = os.path.join(cnv_out,"manta")
		manta_out = os.path.join(temp_dir,"manta")
	
		if not os.path.isdir(manta_out):
			os.makedirs(manta_out)
	
		# NOTE - THIS SHOULD BE CREATED FROM PANEL-CSV IN THE FUTURE
		# Argh, can only do region by chromosome (instead of coordinates) due to stupid filename being too long
	
		manta_target = ""
		manta_target = manta_target + "--region chr1:45794900-45806000 " #MUTYH
		manta_target = manta_target + "--region chr2:47596000-47710100 " #EPCAM-MSH2
		manta_target = manta_target + "--region chr2:48010000-48034000 " #MSH6
		manta_target = manta_target + "--region chr2:215593300-215674400 " #BARD1
		manta_target = manta_target + "--region chr3:37035000-37092200 " #MLH1
		manta_target = manta_target + "--region chr5:112090500-112180000 " # APC
		manta_target = manta_target + "--region chr7:6013000-6048700 " # PMS2
		manta_target = manta_target + "--region chr8:90947700-90996900 " # NBN
		manta_target = manta_target + "--region chr9:97863900-98011700 " # FANCC
		manta_target = manta_target + "--region chr10:89624000-89725200 " # PTEN
		manta_target = manta_target + "--region chr11:108098300-108236300 " #ATM
		manta_target = manta_target + "--region chr13:32890500-32973000 " #BRCA2 
		manta_target = manta_target + "--region chr16:23614700-23652500 " #PALB2
		manta_target = manta_target + "--region chr16:68771200-68867500 " #CDH1
		manta_target = manta_target + "--region chr17:7572900-7580000 " #TP53
		manta_target = manta_target + "--region chr17:33427900-33446700 " #RAD51D
		manta_target = manta_target + "--region chr17:41197600-41276200 " #BRCA1
		manta_target = manta_target + "--region chr17:56769900-56811700 " #RAD51C
		manta_target = manta_target + "--region chr17:59760600-59939000 " #BRIP1
		manta_target = manta_target + "--region chr19:1206900-1226600 " # STK11
		manta_target = manta_target + "--region chr22:29083800-29130800" #CHEK2 # No space for the last gene on the list.
		
		"""
		manta_target = ""
		manta_target = manta_target + "--region chr1 " #MUTYH
		manta_target = manta_target + "--region chr2 " #EPCAM,MSH2,MSH6
		manta_target = manta_target + "--region chr3 " #MLH1
		manta_target = manta_target + "--region chr5 " # APC
		manta_target = manta_target + "--region chr7 " # PMS2
		manta_target = manta_target + "--region chr10 " # PTEN
		manta_target = manta_target + "--region chr11 " #ATM
		manta_target = manta_target + "--region chr13 " #BRCA2 
		manta_target = manta_target + "--region chr16 " #PALB2
		manta_target = manta_target + "--region chr17 " #BRCA1, TP53
		manta_target = manta_target + "--region chr19 " # STK11
		manta_target = manta_target + "--region chr22" #CHEK2 # No space for the last gene on the list.
		"""
	
		# Manta - Builds a python workflow.py
	
		#print (" ".join([mantaCMD, bam_list_file, manta_out, genome, '"' + manta_target + '"', '"' + bam_suffix + '"']))
	
		# The inserted -- stops python argparse from looking for additional optional arguments, and therefore allows for positional arguments to start with a -
		subprocess.call(" ".join([mantaCMD, "--",bam_list_file, manta_out, genome, '"'+manta_target+'"',bam_suffix]), shell=True)
	
		# Clean up intermediate / old files
		manta_workspace_folders = glob.glob(os.path.join(manta_out, "*", "workspace"))
		for workspace_dir in manta_workspace_folders:
			#workspace_dir = os.path.dirname(workspace)
			print ("Removing temporay workspace " + workspace_dir)
			shutil.rmtree(workspace_dir)
	
		temp_ba_files  = glob.glob(os.path.join(temp_dir, "*.ba*"))
	
		for ba_file in temp_ba_files:
			print ("Removing temporary BAM/BAI " + ba_file)
			os.remove(ba_file)
	
		# Copy all manta output vcf to run's cnv folder
		manta_call_files = glob.glob(os.path.join(manta_out, "*manta.cnv.vcf"))
		manta_final_out = os.path.join(cnv_out,"manta")
	        
	        if not os.path.isdir(manta_final_out):
	                os.makedirs(manta_final_out)
	
		for call_file in manta_call_files:
			print ("Copying " + call_file + " to Alignment Folder and cleaning temp folder")
			shutil.copy(call_file, manta_final_out)
			os.remove(call_file)
	
#Defining internal functions:
def log_write(log_file, msg):
	log_handle=open(log_file, 'a')
	log_handle.write (msg)
	log_handle.close()

main()
	
	
