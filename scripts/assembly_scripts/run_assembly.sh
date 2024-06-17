NAME=$1
HiFi=$2
#HC1=$3
#HC2=$4

echo "Processing ${ASSEMBLY}"
#### Edit ploidy here if needed 
#l_value="3"

echo $HiFi

###Activate conda env for packages
#conda activate rice_pang		#activate conda env before running the script

###run KMC:
mkdir ${NAME} ${NAME}_kmc
srun -c 80 kmc -k21 -m90 -t80 ${HiFi} ${NAME} ${NAME}_kmc
kmc_tools -t20 transform ${NAME} histogram ${NAME}_hist.txt

### run hifiasm
mkdir ${NAME}_hifiasm
srun -c 200 hifiasm -t 200 -o ${NAME}_hifiasm -l0 --write-ec ${HiFi} 2> ./${NAME}_hifiasm.stdout
mv ${NAME}_hifiasm.* ./${NAME}_hifiasm

### Fix output gfas to fas
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.bp.p_ctg.gfa > ./${NAME}_hifiasm/${NAME}_hifiasm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.bp.a_ctg.gfa > ./${NAME}_hifiasm/${NAME}_hifiasm.bp.a_ctg.fa
cat ./${NAME}_hifiasm/${NAME}_hifiasm.bp.p_ctg.fa ./${NAME}_hifiasm/${NAME}_hifiasm.bp.a_ctg.fa > ./${NAME}_hifiasm/${NAME}_hifiasm.combined_ctg.fa

### Map reads to prepare for polishing primary and combined haplotype assembly. Hash out if only one is needed.
srun -c 128 minimap2 -x map-pb -I 128G -a -t 128 -c --secondary=no ./${NAME}_hifiasm/${NAME}_hifiasm.combined_ctg.fa ${HiFi} > ./${NAME}_hifiasm/${NAME}_HiFi_combined.sam

### run racon polishing
mkdir ${NAME}_racon
srun -c 100 racon -t 100 ${HiFi} ./${NAME}_hifiasm/${NAME}_HiFi_combined.sam ./${NAME}_hifiasm/${NAME}_hifiasm.combined_ctg.fa > ./${NAME}_racon/${NAME}_hifiasm.combined_racon.fa


### Now assembly assessment
ASSEMBLY="${NAME}_hifiasm.combined_racon"
ASSEMBLY_PATH="./${NAME}_racon"
BUSCO="poales_odb10"

### Run assembly_stats and BUSCO
sh ./scripts/assembly_scripts/assembly_assessment_1.sh ${ASSEMBLY} ${BUSCO} ${ASSEMBLY_PATH}

### Run Merqury
## May need to edit range for plot
sh ./scripts/assembly_scripts/assembly_assessment_2.sh ${ASSEMBLY} ${HiFi} ${ASSEMBLY_PATH}


#### Map reads and create coverage plot of assembly. 
### May need to edit plotting length requirements for in assembly_assessment_3.sh if outside of range
sh ./scripts/assembly_scripts/assembly_assessment_3.sh ${ASSEMBLY} ${HiFi} ${ASSEMBLY_PATH}

### Files housekeeping
mv ${NAME}* ${NAME}/
mv log* ${NAME}/
mv Assembly_Assessment_${NAME}* ${NAME}

echo "####################################################################"
echo "Pipeline run done"
date
echo "####################################################################"
