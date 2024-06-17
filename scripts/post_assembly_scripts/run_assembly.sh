NAME=$1
HiFi=$2
HC1=$3
HC2=$4

echo "Processing ${ASSEMBLY}"
#### Edit ploidy here if needed 
l_value="3"

echo $HC1
echo $HC2
echo $HiFi

mkdir ${NAME}_hifiasm

srun -c 200 /data2/castri/hifiasm-0.19.6/hifiasm -t 200 -o ${NAME}_hifiasm -l${l_value} --write-ec --h1 ${HC1} --h2 ${HC2} ${HiFi} 2> ./${NAME}_hifiasm.stdout
mv ${NAME}_hifiasm.* ./${NAME}_hifiasm

## Fix output gfas to fas
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.gfa > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap1.p_ctg.gfa  > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap2.p_ctg.gfa  > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap2.p_ctg.fa
cat ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap1.p_ctg.fa ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap2.p_ctg.fa > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap12.p_ctg.fa

## Map reads to prepare for polishing primary and combined haplotype assembly. Hash out if only one is needed.
srun -c 128 minimap2 -I 128G -ax map-hifi -t 128 -c --secondary=no ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap12.p_ctg.fa ${HiFi} > ./${NAME}_hifiasm/${NAME}_HiFi2_HAP12.sam
#srun -c 128 minimap2 -I 128G -ax map-hifi -t 128 -c --secondary=no ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.fa ${HiFi} > ./${NAME}_hifiasm/${NAME}_HiFi2_primary.sam

### run racon polishing
mkdir ${NAME}_racon
srun -c 100 /data2/castri/racon/build/bin/racon -t 100 ${HiFi} ./${NAME}_hifiasm/${NAME}_HiFi2_HAP12.sam ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap12.p_ctg.fa > ./${NAME}_racon/${NAME}_hifiasm.HAP12_racon.fa
#srun -c 100 /data2/castri/racon/build/bin/racon -t 100 ${HIFI} ./${NAME}_hifiasm/${NAME}_HiFi2_primary.sam ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.fa > ./${NAME}_racon/${NAME}_hifiasm.primary_racon.fa

### Now assembly assessment

ASSEMBLY="${NAME}_hifiasm.HAP12_racon"
ASSEMBLY_PATH="./${NAME}_racon/"
BUSCO="eudicots_odb10"

#### Run assembly_stats and BUSCO
source activate busco_env
sh assembly_assessment_1.sh ${ASSEMBLY} ${BUSCO} ${ASSEMBLY_PATH}
source deactivate

#### Run Merqury
### May need to edit range for plot
source activate r_env
sh assembly_assessment_2.sh ${ASSEMBLY} ${HiFi} ${ASSEMBLY_PATH}
source deactivate

#### Map reads and create coverage plot of assembly. 
#### May need to edit plotting length requirements for in assembly_assessment_3.sh if outside of range
sh assembly_assessment_3.sh ${ASSEMBLY} ${HiFi} ${ASSEMBLY_PATH}