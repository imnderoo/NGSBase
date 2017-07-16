#/bin/bash

run_folder=$1
target_bed="/media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed"

python /media/sf_resources/tools/create_report_v3/docker_pipeline.py $run_folder/Data/Intensities/BaseCalls/ $target_bed $run_folder/Data/Intensities/BaseCalls/Alignment

python /media/sf_resources/tools/create_report_v3/create_report.py $run_folder 20 20 /media/sf_PLMGenetics/NGS_Analysis_Output/ /media/sf_resources/ --casefolder "/media/sf_CaseData/Molecular/Analytical Data/"

#python ./create_report.py $run_folder 20 20 /media/sf_PLMGenetics/NGS_Analysis_Output/ /media/sf_resources/

