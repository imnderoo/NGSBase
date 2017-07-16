#/bin/bash

#Val 2
#python ./docker_pipeline_20170111b.py /media/sf_PLMGenetics/NextSeqOutput/170303_NB501823_0007_AHJL3TAFXX/Data/Intensities/BaseCalls/ /media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed /media/sf_PLMGenetics/NextSeqOutput/170303_NB501823_0007_AHJL3TAFXX/Data/Intensities/BaseCalls/Alignment
#python ./create_report.py /media/sf_PLMGenetics/NextSeqOutput/170303_NB501823_0007_AHJL3TAFXX/ 20 20 /media/sf_PLMGenetics/Validation/NextSeq/ /media/sf_resources/

#MiSeq Reanalysis - For Variants Comparison
#python ./docker_pipeline_20170111b.py /media/sf_PLMGenetics/MiSeqOutput/170224_M03448_0117_000000000-B338M/Data/Intensities/BaseCalls/ /media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed /media/sf_PLMGenetics/MiSeqOutput/170224_M03448_0117_000000000-B338M/Data/Intensities/BaseCalls/Alignment

#python ./create_report.py /media/sf_PLMGenetics/MiSeqOutput/170224_M03448_0117_000000000-B338M/ 20 20 /media/sf_PLMGenetics/Validation/NextSeq/ /media/sf_resources/

#python ./docker_pipeline_20170111b.py /media/sf_PLMGenetics/MiSeqOutput/161216_M03063_0115_000000000-ATG0P/Data/Intensities/BaseCalls/ /media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed /media/sf_PLMGenetics/MiSeqOutput/161216_M03063_0115_000000000-ATG0P/Data/Intensities/BaseCalls/Alignment

#python ./create_report.py /media/sf_PLMGenetics/MiSeqOutput/161216_M03063_0115_000000000-ATG0P/ 20 20 /media/sf_PLMGenetics/Validation/NextSeq/ /media/sf_resources/

#python ./docker_pipeline_20170111b.py /media/sf_PLMGenetics/MiSeqOutput/160902_M03063_0096_000000000-AT17G/Data/Intensities/BaseCalls/ /media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed /media/sf_PLMGenetics/MiSeqOutput/160902_M03063_0096_000000000-AT17G/Data/Intensities/BaseCalls/Alignment

#python ./create_report.py /media/sf_PLMGenetics/MiSeqOutput/160902_M03063_0096_000000000-AT17G/ 20 20 /media/sf_PLMGenetics/Validation/NextSeq/ /media/sf_resources/

run_folder=$1
target_bed="/media/sf_resources/createReport_DBs/analysis_type/msh21_inherited/msh21_inherited.bed"

#python ./docker_pipeline.py $run_folder/Data/Intensities/BaseCalls/ $target_bed $run_folder/Data/Intensities/BaseCalls/Alignment

python /media/sf_resources/tools/create_report_v3/create_report.py $run_folder 20 20 /media/sf_PLMGenetics/NGS_Analysis_Output/ /media/sf_resources/ --casefolder "/media/sf_CaseData/Molecular/Analytical Data/"

