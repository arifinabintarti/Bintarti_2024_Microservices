Raw Sequence Analysis of Bacterial amoA Gene (AOB) of Wheat Field Soil and Rhizosphere Samples from the Rain-out Shelter and DOK Experiment (MICCROSERVICES)

# Analysis of amoA gene of AOB Illumina MiSeq Data

The sequencing library was constructed using two-step PCR to enable for primer barcoding. Sequencing was conducted at Sequencing Facility of Genoscreen in Lille, France using V2 flow cell with 2 x 250 bp.

All raw sequence analyses were conducted using the INRAE server that can be accessed by ssh connection and INRAE authenticators: 
$>ssh <inrae-username>@djs-agro-emfeed.inra.local

All data are stored in the "microservices" directory under "Project" directory in the server (/home/afbintarti/Projects/microservices). Raw sequence data are stored in the INRAE server:
/home/afbintarti/Projects/microservices/03042023_AOA_AOB_rawSeq/AOB

The pipeline used for the analysis is AMOA-SEQ pipeline developed by Lee et  al. 2023 (it has not published yet), a research team at the University of Lyon, France. The pipeline uses DADA2 for sequence quality filtering to generate amplicon sequence variants (ASVs) (https://github.com/miasungeunlee/AMOA-SEQ).

Because the quality of the AOB raw sequece data is low (the average length of ~89 bp) merging of paired reads was not possible. Thus, we used GAP pipeline with the AOB sequence data.

Follow the tutorial and installation of the AMOA-SEQ tool in: https://github.com/miasungeunlee/AMOA-SEQ.




