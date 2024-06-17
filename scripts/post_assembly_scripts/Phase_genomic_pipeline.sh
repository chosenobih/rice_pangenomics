### Map HiC for Phase Genomics, etc (Not Arima)
### Require bwa, haphic environment and samtools

### After creating sample info list, execute with: sh run_assembly <sample.list>
### See below for sample info list details

SAMPLE_LIST=$1
if [ -z "${SAMPLE_LIST}" ]
then
    echo "This script needs the SAMPLE lists as the first parameter"
    echo "the SAMPLE info list is a CSV with 6 columns"
    echo "the 1st column is sample name which will be the prefix for file names"
    echo "the 2nd column is HiFi reads including path"
    echo "the 3rd column is Hi-C read pair 1 including path"
    echo "the 4th column is Hi-C read pair 2 including path"
    echo "Expected ploidy of sample"
    echo "Expected chromosome number"
    exit
else
    echo "Processing SAMPLEs from ${SAMPLE_LIST}"
fi
 
#---- PROCESSING
 
for SAMPLEINFO in $(cat $SAMPLE_LIST);
 
do
  NAME=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $1}')
  HiFi=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $2}')
  HC1=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $3}')
  HC2=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $4}')
  Ploidy=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $5}')
  Chr=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $6}')
  
ASSEMBLY="../${NAME}_racon/${NAME}_hifiasm.HAP12_racon.fa"
OUTPUT_NAME="${NAME}_Phase"

echo $ASSEMBLY
echo $HC1
echo $HC2

bwa index ./${NAME}_racon/${NAME}_hifiasm.HAP12_racon.fa
srun -c 60 bwa mem -SP -t 60 ${ASSEMBLY} ${HC1} ${HC2} | samtools view -bhS - > ./${NAME}_racon/${OUTPUT_NAME}.bam

HIC_BAM="./${NAME}_racon/${OUTPUT_NAME}.bam"

### Run HapHiC pipeline
source activate haphic
mkdir ${NAME}_HapHIC
cd ${NAME}_HapHIC
srun -c 100 /data2/castri/HapHiC/haphic pipeline ${ASSEMBLY} ${HIC_BAM} ${CHR} --correct_nrounds 2 --remove_allelic_links ${PLOIDY} --threads 50 --prefix ${NAME}_Haphic

### Create files for Juicebox
cd ./04.build 
bash juicebox.sh
cd ../..
source deactivate
done
