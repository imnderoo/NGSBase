#/bin/bash

run_folder=$1
resource_folder="/mnt/resources/"
target_bed=$resource_folder/createReport_DBs/analysis_type/acmg59/acmg59.bed
out_folder="/mnt/analysis_output"

python $resource_folder/tools/create_report_v3/docker_pipeline.py $run_folder/Data/Intensities/BaseCalls/ $target_bed $run_folder/Data/Intensities/BaseCalls/Alignment

python $resource_folder/tools/create_report_v3/create_report.py $run_folder 20 20 $out_folder $resource_folder
