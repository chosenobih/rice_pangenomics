#### Genome Assembly README

The genome assembly pipeline workflow of scripts goes as the following:

1) run_assembly.sh 
2) Arima_pipeline.sh or Phase_genomic_pipeline.sh
3) add_organelles.sh
4) contamination_screen.sh
5) EDTA_tRNA.sh
6) braker.sh

1) To run run-assembly.sh: sh run_assembly.sh <sample info list>
   The input the sample info list is a CSV with 6 columns
   -the 1st column is sample name which will be the prefix for file names"
   -the 2nd column is HiFi reads including path"
   -the 3rd column is Hi-C read pair 1 including path"
   -the 4th column is Hi-C read pair 2 including path"
   -the 5th columns is the expected ploidy of sample"
   -the 6th column is the expected chromosome number"
   
   
If there are multiple samples, put them as new rows in the list and it will run them one by
one. 

This script first creates a hifiasm assembly using hifi and hic reads. In the script the 
default ploidy is set to 2+, if the ploidy is 1 this value should be edited in the script.
Next the output files are changed to .fa instead of gfa. 

Then the raw hifi reads are used to polish the genome using racon. The script is currently
set up to utilize the diploid assembly instead of the separate haploid ones. If the primary
assembly is desired, unhash out that line and hash out the diploid line. 

Lastly, the script runs the assembly assessment workflow on the polished assembly. The 
assembly assessment workflow consists of assembly_stats statistics, busco analysis,
Merqury analysis and plot, coverage plotting with wgsplotter. There are three 
assembly_assessment_1-3.sh scripts that split up the analyses. There is also a run_assembly_assessment.sh
script that is not part of the initial assembly pipeline, but can be used in later stages whenever 
an assembly should be assessed without the prior hifiasm and racon. The three scripts for
the assessment are to make it easy to run only parts of the assessment and all that is required 
is to hash out the undesired ones in the run_assembly_assessment.sh script if this is
wanted. 

To run the assembly_assessment pipeline at later stage of the assembly process do the 
following:

sh run_assembly_assessment.sh <sample_info>
The path and name of new assembly need to be updated in the run_assembly_assessment.sh script for later
<assembly name without the .fa the end!> <path to assembly, no / at the end)
The default BUSCO database in both the run_assembly.sh and run_assembly_assessment scripts is
eudicots_odb10. If another database is needed, this needs to be changed in the scripts.

The following dependencies are needed for this script:
hifiasm (currently it is calling the newest version without a bug in my folder. Change if needed)
minimap2
racon (Again, it is calling the installation from my directory. Change if needed)
assembly_stats
conda busco installation (mine is called busco_env)
conda R environment with Merqury dependencies installed (script calls r_env)
meryl
Merqury
bwa
samtools
gatk docker container
wgsplotter (called from my directory. This requires java. module load java/21 is included in script)

The Merqury plot and wgsplot can be finicky and the axes need to be readjusted. Wgsplot
will sometimes fail if the contig size to plot is too small or big. 

2) The HiC mapping step is either Arima_mapping_pipeline.sh or Phase_genomics_pipeline.sh
depending on the data type. 

To run:
sh Phase_genomics_pipeline.sh <sample.list>
sh Arima_pipeline.sh <sample.list>

These scripts map the HiC reads to the assembly depending on the type of kit used for preparing the reads.
The Phase genomics uses bwa to map HiC reads while the Arima reads are mapped using the Arima
mapping pipeline that they provide. Included in this is read trimming and a number of processing steps.
The last step in both scripts is to scaffold into pseudomolecules using HapHic and the produced bam files.
Once the new assembly is created, the input files for Juicebox manual inspection/editing are
created. The next step is to upload the .assembly and .hic files into Juicebox for manual editing.

The following dependencies are needed for these scripts:
bwa
samtools
java/21
perl
picard (script calls from my installation)
HapHiC conda environment, mine is called haphic

3) The script for adding organelles uses oatk to identify and build organelles. The script
first extracts the primary assembly and unanchored scaffolds from the current assembly. Oatk
is run on the unanchored scaffolds and then added to the primary assembly. The unanchored scaffolds
are then mapped to the primary assembly and organelles to determine which belong to the orgnanelles.
The organelle scaffold are then removed and the final assembly is the primary assembly, organelles,
unplaced scaffolds that don't map to organelles.  This script is still in development and may require
tweaking. Oatk is also a new program and I have had some weird results where it named the
mito and pltd contigs the same name which created issues, while it worked perfectly other times.
The database has been run on the angiosperm mitochondria and chloroplast databases which were
downloaded from oatk's github. 

To run: sh add_organelles.sh <sample.list>

Dependencies:
oatk (and mito and pltd databases)
minimap2
samtools
seqtk

4) Contamination screen is preformed with fcs-gx. It requires the taxon ID for the sample which
can be found on the NBCI webpages. It also requires the ncbi database is installed, the script currently points 
to my installation. The program fails with other python versions, so run with python37. It also fails when run 
with srun unfortunately so it can take a bit.

To run: sh contamination_screen.sh <assembly_file> <output_directory_name> <NBCI ID>

Dependencies:
python37
fcs-gx (including its gxdb database)

5) The EDTA and tRNAscan script runs both but one can be hashed out or the other if only one
is needed. A new csv file is needed for this script which includes the samplename, assembly including path, 
and species which should be one of the following [Rice,Maize,others]. 

To run: sh EDTA_tRNA.sh <sample.list>

Dependencies:
EDTA
tRNAscan-SE

6) The braker.sh script is the least automated in the pipeline as I've only ran braker once
and the input data types are expected to be different for other samples. I thought it would be
annoying to make another input info data sheet and the script is likely going to need to be
edited anyway for instance, if there is no Isoseq data, if there is only Isoseq etc. Just hash out
the variables and script lines that are not applicable to the sample. All input information is edited
in the beginning of file. 

To run: sh braker.sh

Dependencies:
ULTRA
histat2
minimap2
samtools
braker singularity container (regular container, for Isoseq container (braker_lr.sif), or both)
tsebra
gtf2aa.pl from Augustus (if you want to extract proteins from gtf. gffread does not parse correctly. 
   Available at https://fossies.org/dox/Augustus-3.5.0/gtf2aa_8pl_source.html)


