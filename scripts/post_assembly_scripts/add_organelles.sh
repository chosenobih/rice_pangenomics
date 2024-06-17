### Identify and extract organelles from assembly then incorporate into final assembly 
### Requires: oatk and appropriate organelle database (in working directory),
### FASTA of pseudomolecules (CHR_ASSEMBLY), full assembly file, and unanchored scaffold FASTA file. Needs to be done manually once CHRs are determined
### To run sh add_organelles.sh <sample.list>
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
  

######## Edit info below #########
OUTPUT="${NAME}"
MITO="angiosperm_mito.fam"
CP="angiosperm_pltd.fam"


ASSEMBLY="./${NAME}_HapHIC/04.build/${NAME}.FINAL.fa"

samtools faidx ${ASSEMBLY}

sort -nk 2 ${ASSEMBLY}.fai | cut -f 1 | head -n ${Chr} > primary_scaffolds.list
scaffold_count=`wc -l ${ASSEMBLY}.fai |  awk '{ print $1 }'`
unanchored_count=`echo $((${scaffold_count}-${Chr}))`
sort -nk 2 ${ASSEMBLY}.fai | cut -f 1 | tail -n ${unanchored_count} > unanchored_scaffolds.list

seqtk subseq ${ASSEMBLY} primary_scaffolds.list > ./${NAME}_HapHIC/${NAME}.FINAL.primary.fa
seqtk subseq ${ASSEMBLY} unanchored_scaffolds.list > ./${NAME}_HapHIC/${NAME}.FINAL.unanchored.fa

CHR_ASSEMBLY="./${NAME}_HapHIC/${NAME}.FINAL.primary.fa"
UNANCHORED="Cereus_jamacaru_final_unanchored.fa"

srun -c 50 oatk -k 1001 -c 30 -t 50 -m ${MITO} -p ${CP} -o ${OUTPUT} ${HiFi}

######  If there is an issue with the pipeline, these are subscripts that can be used to initate the later stages
#hmm_annotation -t 8 -o ${OUTPUT}.annot_mito.txt ${MITO} ${OUTPUT}.utg.final.gfa
#hmm_annotation -t 8 -o ${OUTPUT}.annot_pltd.txt ${CP} ${OUTPUT}utg.final.gfa
#pathfinder -m ${OUTPUT}.annot_mito.txt -p ${OUTPUT}.annot_pltd.txt -o ${OUTPUT} ${OUTPUT}.utg.final.gfa

cat ${CHR_ASSEMBLY} ${OUTPUT}.mito.ctg.fasta  ${OUTPUT}.pltd.ctg.fasta > ${OUTPUT}_primary_organelles.fasta

## Map back unanchored scaffolds to organelles to extract mt and pltd scaffolds
srun -c 128 minimap2 -I 128G -ax map-hifi -t 128 -c --secondary=no ${OUTPUT}_primary_organelles.fasta ${UNANCHORED} > ${OUTPUT}_unanchored2primary_organelle.sam

### Get names of organelles from oatk
grep ">" ${OUTPUT}.mito.ctg.fasta | sed 's/>//g' | cut -f 1 > mito.list
grep ">" ${OUTPUT}.pltd.ctg.fasta | sed 's/>//g' | cut -f 1 > pltd.list
cat mito.list chloroplast.list > organelles.list

### Prepare bam file
samtools sort ${OUTPUT}_unanchored2primary_organelle.sam -o ${OUTPUT}_unanchored2primary_organelle.bam
samtools index ${OUTPUT}_unanchored2primary_organelle.bam

### Extract scaffolds mapping to organelles
samtools view -q 20 ${OUTPUT}_unanchored2primary_organelle.bam organelles.list > ${OUTPUT}_organelles.sam
### Get list of scaffolds mapped to organelles
cut -f 1 ${OUTPUT}_organelles.sam > organelle_scaffolds.list
### Get list of all initial assembly scaffolds
grep ">" ${ASSEMBLY} | sed 's/>//g' | cut -f 1 | sort > all_scaffolds.list
#### Remove organelle scaffolds from all assembly scaffolds so they are not in assembly twice
comm -3 all_scaffolds.list organelle_scaffolds.list > scaffolds_no_organelles.list

### Create assembly FASTA without organelles or the scaffolds aligned to them
seqtk subseq ${ASSEMBLY} scaffolds_no_organelles.list > ${ASSEMBLY}_no_organelles.fa

### Combine assembly with no organelles with the oatk FASTA files of organelles for final assembly
### Stupid name sorry, rename accordingly
cat ${ASSEMBLY} ${OUTPUT}.mito.ctg.fasta ${OUTPUT}.pltd.ctg.fasta > ${OUTPUT}_final_with_organelles.fa
