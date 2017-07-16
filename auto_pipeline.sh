#!/bin/bash

for run in $(find /media/sf_PLMGenetics/MiSeqOutput/ -maxdepth 1 -newerct "mar 31 2017" -name "*")
do
	if [ -e $run/SampleSheet.csv ]
	then
		if grep -q "Nextera Rapid Capture Enrichment" $run/SampleSheet.csv || grep -q "TruSight Enrichment" $run/SampleSheet.csv
		then 	
			if [ -e $run/CompletedJobInfo.xml ] && [ -e $run/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz ]
			then
				if [ ! -e $run/NGS_Analysis_Complete.txt ] && [ ! -e $run/NGS_Analysis_Running.txt ]
				then
					echo "$run not completed. FASTQ detected."
					echo "Auto-starting Analysis: $run"
					touch $run/NGS_Analysis_Running.txt
		
					#find $run/Data/Intensities/BaseCalls/ -maxdepth 1 -name "*fastq.gz" > $run/auto_find_fastq.txt

					/media/sf_resources/tools/create_report_v3/docker_pipeline_analysis.sh $run 2>&1 | tee $run/analysis_log.txt
					
					rm -f $run/NGS_Analysis_Running.txt
					touch $run/NGS_Analysis_Complete.txt
					
				fi
			fi
		fi
	fi
done
