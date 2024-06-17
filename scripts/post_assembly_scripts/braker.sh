### Annotation workflow with IsoSeq, short RNA and protein database using Braker3
### Requires braker3 singularity container and braker2_lr singularity container, protein database (orthodb etc)
### also hisat, minimap2, TSEBRA, ultra, EDTA, tRNAScan-SE are required 

### ASSEMBLY is the prefix here. Fasta should have ".fa" extension




ASSEMBLY_PREFIX=""
ASSEMBLY_PATH=""
### ISOSEQ READS FASTA
ISOSEQ=""
NAME=""
SHORT_1=""
SHORT_2=""
PROT_DB=""
CHR_PREFIX=""

### Make tandem repeats
srun -c 50 ultra -o ${ASSEMBLY_PATH}/${ASSEMBLY_PREFIX}.bed -t 50 --mask ${ASSEMBLY}/${ASSEMBLY_PREFIX}_masked.fa ${ASSEMBLY}/${ASSEMBLY_PREFIX}.fa

### Map IsoSeq
minimap2 -ax splice:hq -uf ${ASSEMBLY}/${ASSEMBLY}_masked.fa ${ISOSEQ} > IsoSeq2${ASSEMBLY}.sam    # Final PacBio Iso-seq or traditional cDNA

### Map short read RNA
srun -c 50 BuildDatabase --threads 50 -name ${NAME} ${ASSEMBLY}/${ASSEMBLY_PREFIX}.fa

hisat2-build ${ASSEMBLY}/${ASSEMBLY_PREFIX}_masked.fa ${ASSEMBLY}/${NAME}_masked
hisat2 --dta -p 20 --max-intronlen 100000 -x ${ASSEMBLY}/${NAME}_masked -1 ${SHORT_1} -2 ${SHORT_2}  -S ${NAME}_rna.sam
samtools sort ${NAME}_rna.sam -o ${NAME}_rna.bam

### Braker pipeline with long and short RNA and protein db

srun -c 100 singularity.AGI exec -B ${PWD}:${PWD} braker3_lr.sif braker.pl --genome=${ASSEMBLY}/${ASSEMBLY_PREFIX}_masked.fa --bam=IsoSeq2${ASSEMBLY}.bam --threads=100
mv braker braker_longRNA
srun -c 100 singularity.AGI exec -B ${PWD}:${PWD} braker3.sif braker.pl --genome=${ASSEMBLY}/${ASSEMBLY_PREFIX}_masked.fa --prot_seq=${PROT_DB} --bam=${NAME}_rna.bam --threads=48
mv braker braker_shortRNA_protein

### Combine gtfs from braker runs
/data2/castri/TSEBRA/bin/tsebra.py -e ./braker_shortRNA_protein/hintsfile.gff,./braker_longRNA/hintsfile.gff -g braker_shortRNA_protein/augustus.hints.gtf,braker_longRNA/braker.gtf -c /data2/castri/TSEBRA/config/default.cfg \
-o RNA_short_long_protein_filtered.gtf --filter_single_exon_genes

### Add prefix to chromosomes
/data2/castri/TSEBRA/bin/rename_gtf.py --gtf RNA_short_long_protein_filtered.gtf --prefix ${PREFIX}

### Create proteins fasta from gtf to check with BUSCO for example
gtf2aa.pl $/${ASSEMBLY}.fa RNA_short_long_protein_filtered_renamed.gtf ${NAME}_proteins.fa
