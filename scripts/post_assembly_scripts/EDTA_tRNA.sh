### Run EDTA and tRNAscan

SAMPLE_LIST=$1

if [ -z "${SAMPLE_LIST}" ]
then
    echo "This script needs the SAMPLE lists as the first parameter"
    echo "the SAMPLE info list is a CSV with 3 columns"
    echo "the 1st column is sample name which will be the prefix for file names"
    echo "the 2nd column is genome assembly including path"
    echo "the 3rd column is species to run the curated library for [Rice,Maize,others]"

    exit
else
    echo "Processing SAMPLEs from ${SAMPLE_LIST}"
fi
 
#---- PROCESSING
 
for SAMPLEINFO in $(cat $SAMPLE_LIST);
 
do
  NAME=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $1}')
  ASSEMBLY=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $2}')
  SPECIES=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $3}')

echo "Processing ${NAME}"

### Run EDTA
srun -c 50 perl /data2/castri/EDTA/EDTA.pl --genome ${ASSEMBLY}.fa --threads 50 --species ${SPECIES}
### Run tRNAscan-SE
srun -c 100 tRNAscan-SE --thread 100 -j ${NAME}_trna.gff -a ${NAME}_trna.fa -o ${NAME}_trna_results -m ${NAME}_trna_stats.txt ${ASSEMBLY}.fa 
