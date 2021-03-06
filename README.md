# STMAP
Installation Instructions
================================
  
To install Metagenome Analysis Pipeline, you must have a minimum of 8 GiB free disk space and minimum of 16 GiB free RAM to test run. 

To provide an easier way to install, we provide a miniconda based installer.
Installation also requires **pre-instaled** ``git``, ``gcc``, ``cpp`` and ``zlib1g-dev``.
  
    git clone https://github.com/computational-genomics-lab/STMAP.git
    cd STMAP
    chmod 755 INSTALL.sh
    ./INSTALL.sh

    
**Post Installation Instructions**

After successful installation, close the current terminal. 
In a new terminal. source the bashrc file:  ``source ~/.bashrc``
One note of caution here is if your .bashrc file is already populated and is tuned to run some of our other pipelines, then clean it 
up before installing STMAP. It may interfere with the environment variables.

All the third party tools installed using conda are available at $HOME/STMAP/ [default location]
or the user specified location during the installation process.

The script to run Metagenome Analysis Pipeline is metagenome.py is available inside the STMAP folder, that you cloned from github.

**Installation of enrichM database**

Please follow the instruction provided at https://github.com/geronimp/enrichM#setup to install and set path to enrichM database.

**After installing all pre-rquisites for enrichm**
After installation of enrichm and its other requisite packages, the next thing to be installed is the database using the following command

    enrichm data
    
  **enrichM database Installation Issue**
  
  [However, it is noted that many times, this command does not run well. In case you 
  encounter an error as below then]




  File "$HOME/STMAP/envs/enrichm/bin/enrichm", line 357, in <module>
    r.main(args, sys.argv)
  File "$HOME/STMAP/envs/enrichm/lib/python3.7/site-packages/enrichm/run.py", line 290, in main
    d.do(args.uninstall)
  File "$HOME/STMAP/envs/enrichm/lib/python3.7/site-packages/enrichm/data.py", line 120, in do
    version_remote = urllib.request.urlopen(self.ftp + self.VERSION).readline().strip().decode("utf-8")
AttributeError: module 'urllib' has no attribute 'request'

  Go to files:
  1. $HOME/STMAP/envs/enrichm/lib/python3.7/site-packages/enrichm/data.py
  2. $HOME/STMAP/envs/enrichm/lib/python3.7/site-packages/enrichm/run.py
  
  add 'import urllib.request' (without quotes)
  run
    enrichm data    
   
   again ( At this stage a data of size 5.6 GB is downloaded and may take approximately 15 mins)
  
  

Input Files
===========

Raw Reads
----------------

  The raw ilumina reads in FASTQ format must be placed inside a folder with read permission
  Allowed ``extension`` for FASTQ reads: ``fq`` , ``fq.gz`` , ``fastq``, ``fastq.gz``
 
   For **paired-end** RNAseq reads, sample name must be suffixed with _R1. ``extension`` and _R2. ``extension`` for forward and reverse reads respectively

          *Example*

           sample1_R1.fastq 
           sample1_R2.fastq

           sample2_R1.fastq.gz 
           sample2_R2.fastq.gz

           sample3_R1.fq
           sample3_R2.fq

           sample4_R1.fq.gz 
           sample4_R2.fq.gz
          
   where  ``sample1`` , ``sample2`` , ``sample3`` , ``sample4`` are the sample names
           sample_name must be suffixed with _R1.{extension} and _R2.{extension}

   
The illumina paired end sample data for testing the pipeline can be obtained from
https://figshare.com/articles/dataset/Sample_Data_for_MGAPipe/13340522


    Sample meta data mapping

    samples	     conditions
    SRS043663	   tongue_dorsum
    SRS020336	   buccal_mucosa
    SRS016342	   tongue_dorsum
    SRS016086	   tongue_dorsum
    SRS022532	   buccal_mucosa
    SRS014888	   tongue_dorsum
    SRS017713	   tongue_dorsum
    SRS019327	   tongue_dorsum
    SRS019221	   buccal_mucosa
    SRS013506	   buccal_mucosa
    SRS045049	   buccal_mucosa
    SRS022145	   buccal_mucosa
    SRS019219	   tongue_dorsum
    SRS019329	   buccal_mucosa


Commands
========

**1. Prepare Project**

To design and perform a Metagenome Analysis, a Project need to be prepaired using the ``projectConfig.py`` script. Conda environment must be activated before running the script.

Usage:  projectConfig.py -h

    prepareProject.py <arguments>
    -h       --help             Show this help message and exit

    mandatory arguments         Description
    
    -i      --inputDir          Path to Directory containing raw RNASeq reads, annotation file (gff / gtf), genome and (or) 
                                transcriptome FASTA file. Make sure all the input files are located here.
                                type: string
                                Example: $HOME/STRAP-main/sample_data/mastigocladus

    -e     --emailAddress       Provide your email address
                                type: string
                                Default: Null
                                
    -d     --domain             Organism Domain
                                type: string
                                allowed values: [prokaryote, eukaryote]
                                Default: prokaryote

    
    Optional arguments
    ------------------
    
    -p     --projectName        Name of the Project Directory to be created
                                Must not contain blank spaces and (or) special characters
                                type: string
                                Default: metagenome_project

    -s     --schedulerPort      Scheduler Port Number for luigi
                                type: int
                                Default: 8082
                                
    -o     --symLinkDir         Name of the Symbolic Link Directory to be created
    
    -t     --threads            Number of threads to be used
                                type: int
                                Default = (total threads -1)

    -x     --maxMemory          Maximum allowed memory in GB. 
                                type: int
                                Default = [(available memory in GB) -1)
    

**Run Example**
You can launch a run for projectConfig.py using the following commandline:
    
    mkdir MetagenomeAnalysis
    cd MetagenomeAnalysis
   
    [MetagenomeAnalysis]$ projectConfig.py -i /home/sutripa/MetagenomeAnalysis/data -p metagenome_demo_analysis -d prokaryote -o metagenome_symlink -e tsucheta@gmail.com


Make sure your input directory contains all the data including the metagenome WGS files. Once this commandline is executed, the script is going to iterate over the files in the directory and asks for user input such as whether the read type is pe [Illumina paired end], ont [nanopore] or pac [pacBio] type or if the user wants to exclude the files from data analysis. The workflow takes 2 conditions at a time and generates conditions based on the file pre-fixes. 

If user wants to perform enrichment analysis, then user has to provide conditions associated with each sample.When providing conditions as input make sure to only provide alphanumeric characters. Providing a "-" or "_"  may not work. In case, the file name pre-fixes are in-correct, then one needs to fix it at this time. Also make sure you are not including more than one condition. In that case, the program will exit with error.

If user wants to perform genome resolved metagenomics or taxonomy profiling, user need to provide a single condition name.


Output
    1.  a project folder in the name of ``metagenome_demo_analysis`` will be generated

    2. a configuration folder in the name of config containing 3 files   

      a. metagenome_condition.tsv
      b. metagenome_group.tsv
      b. pe_samples.lst
      c. samples.txt
      
    3. A folder named ``metagenome_symlink`` (provided as parameter to --symLinkDir) will be created which contains the symbolic links to the read files present    in the inputData (/home/sutripa/Documents/scriptome_metaphlan/data) folder

    The ``metagenome_condition.tsv`` file contains the sample names with their associated environmental conditions, which will be used for condition based genome   binning and differential enrichment analysis. Kindly check the generated files and modify if required

Commands to run Metagenome Analysis Workflow
--------------------------------------------
 

metagenome.py ``command`` --help

.. code-block:: none

    Command                      Description   

    1. rawReadsQC                   Raw Reads Quality Assessment using FASTQC tool 

    2. cleanReads                   Process Raw Reads using BBDUK

    3. metagenomeAssembly           Assemble Metagenome using MegaHIT

    4. genomeBinning                Binning individual assembled genome using metabat2
    
    5. binRefinement                Refine Bins using refineM

    6. dRepBins                     de-replicate genome bins using dRep

    7. annotateGenomeBins           Annotation of de-replicated bins using PROKKA

    8. deaGenomeBins                Differential Enrichment Analysis of genome bins (2 different conditions) using enrichM



Prepare Project
===============

To design and perform a Metagenome Analysis experiment, a Project need to be prepaired using the ``configureProject`` command.
The current working directory must contain a template of luigi.cfg file. (https://github.com/computational-genomics-lab/STMAP/blob/main/luigi.cfg)


   Options:
  -h, --help            show this help message and exit


    


  

   **Raw reads quality assessment**

    Before running any of the Metagenome Analysis Workflow commands, a project must be configured using ``configureProject`` command.
    The parent forlder must have the luigi.cfg file, in which the globalparameters are defined.
    Running any of the  Metagenome Analysis Workflow commands without generating the project folder and luigi.cfg file will give rise to     ``luigi.parameter.MissingParameterException``


   **Steps**
   
    1. Run Prepare Projcet with project name ``metagenome_demo_analysis`` as discussed before 
       and inspect the files generated inside config folder

    2. Run rawReadsQC
       [MetagenomeAnalysis]$ metagenome.py rawReadsQC --local-scheduler

      Successful execution of rawReadsQC will generate a folder $PWD/metagenome_demo_analysis/ReadQC/PreQC
      which contains the FASTQC reports of the raw paired-end fastq files


   **Raw samples quality control**
     
    Quality control analysis of the raw samples can be done using command ``preProcessSamples``

    Requirements
       1. Execution of prepareProject.py command 
       2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config``.


    [MetagenomeAnalysis]$ metagenome.py cleanReads <arguments> --local-scheduler
    
    arguments               type      Description
      

     --cleanFastq-kmer-length        int   Kmer length used for finding contaminants
                                      Examle: 13  
                                      Default: 11 

    --cleanFastq-k-trim              str   Trimming protocol to remove bases matching reference
                                      kmers from reads. Choose From['f: dont trim','r: trim to right','l: trim to left] 
                                      Choices: {f, r, l} 

    --cleanFastq-quality-trim       int   Trim read ends to remove bases with quality below trimq.
                                      Performed AFTER looking for kmers.  Values: 
                                          rl  (trim both ends), 
                                          f   (neither end), 
                                          r   (right end only), 
                                          l   (left end only),
                                          w   (sliding window).

                                          Default: f


    --cleanFastq-trim-quality        float  Regions with average quality BELOW this will be trimmed,
                                     if qtrim is set to something other than f.  Can be a 
                                     floating-point number like 7.3  
                                     Default: 6                 

    --cleanFastq-min-length        int   Minimum read length after trimming
                                      Example: 50
                                      Default:50

    --cleanFastq-trim-front          int   Number of bases to be trimmed from the front of the read
                                      Example: 5
                                      Default: 0

    --cleanFastq-trim-tail           int   Number of bases to be trimmed from the end of the read
                                      Example: 5
                                      Default: 0

    --cleanFastq-min-average-quality int   Minimum average quality of reads.
                                      Reads with average quality (after trimming) below 
                                      this will be discarded
                                      Example: 15
                                      Default: 10

    --cleanFastq-mingc              float Minimum GC content threshold
                                      Discard reads with GC content below minGC
                                      Example: 0.1 
                                      Default: 0.0

    --cleanFastq-maxgc              float Maximum GC content  threshold
                                      Discard reads with GC content below minGC
                                      Example: 0.99 
                                      Default: 1.0
    --local-scheduler


**Example Run**

    [MetagenomeAnalysis]$ metagenome.py cleanReads \
                            --cleanFastq-min-average-quality 15 \
                            --cleanFastq-mingc 0.20 \
                            --cleanFastq-maxgc 0.70 \
                            --cleanFastq-quality-trim w  \
                            --local-scheduler

      **Output**
      /path/to/ProjectFolder/metagenome_demo_analysis/ReadQC/CleanedReads/PE-Reads --contains the processed FastQ-reads

A. Genome Resolved Metagenimics
-------------------------------------

**1. Metagenome Assembly using megaHIT**

        Assembly of individual samples can be done using command ``metagenomeAssembly``

        Requirements
        1. Pre execution of ``configureProject`` command 
        2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.



.. code-block:: none   

    [MetagenomeAnalysis]$ metagenome.py  <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the raw illumina reads required or not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

   --local-scheduler




  
  **Metagenome assembly with read quality control analysis** 

     [MetagenomeAnalysis]$ metagenome.py  metagenomeAssembly  \
                                       --pre-process-reads  ``yes``                                    
                                       --local-scheduler



  **Metagenome assembly with out read quality control analysis** 


    [MetagenomeAnalysis]$ metagenome.py  metagenomeAssembly  \
                                       --pre-process-reads  ``yes``                                    
                                       --local-scheduler

     Note: Output folder in the name of MGAssembly will be created at $PWD/metagenome_demo_analysis/MGAssembly
     MGAsembly folder contains the sub-folders by the name of the samples and each sub-folder contains the resultant assembly


**2. Binning of assembled contigs**

Binning of individual assembled samples can be done using command ``genomeBinning``


    Requirements
  
    1. Pre execution of ``configureProject`` command 
    2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
  
  
    [MetagenomeAnalysis]$ metagenome.py  genomeBinning <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the ilumina reads required or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

   --local-scheduler


 
  **Metagenome assembly with read quality control analysis followed by genome binning** 


    [MetagenomeAnalysis]$ metagenome.py  genomeBinning  \
                                       --pre-process-reads  ``yes``                                    
                                       --local-scheduler

    Note: Output folder in the name of binning will be created at ``$PWD/metagenome_demo_analysis/binning``
  ``binning`` folder contains the sub-folders by the name of the samples and each sub-folder contains the resultant genome bins



**3. Refinement of genome bins**

    Refinement of individual genome bins can be done using command ``binRefinement``


    Requirements
  
    1. Pre execution of ``configureProject`` command 
    2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
  
 
    [MetagenomeAnalysis]$ metagenome.py  binRefinement <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the ilumina reads required or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

   --local-scheduler


  **Metagenome assembly with read quality control analysis followed by genome binning and bin refinement** 


    [MetagenomeAnalysis]$ metagenome.py  binRefinement  \
                                       --pre-process-reads  ``yes``                                    
                                       --local-scheduler

    Note: Output folder in the name of binning will be created at ``$PWD/metagenome_demo_analysis/bin_refinement``
  ``bin_refinement``folder contains the sub-folders by the name of the samples and each sub-folder contains the resultant refined genome bins. refineM tool is used      for refinement of genome bins



**4. de-replication of genome bins**

    Condition-based or condition-free de-replication of individual genome bins can be done using command ``dRepBins``
    
    Requirements
    1. Pre execution of ``configureProject`` command 
    2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
    3. Availability of ``metagenome_group.tsv`` file inside the ``config`` folder if user opts for condition_based de-replication
  
    [MetagenomeAnalysis]$ metagenome.py  dRepBins <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the ilumina reads required or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

    --checkM-method        str       choice of workflow to run checkm. choice [lineage_wf, taxonomy_wf]
                                     default: taxonomy_wf

    --dRep-method          str       Mandatory Parameter. choice of de-replication method. choice [condition_based, condition_free]
                    
    --contamination        int       accepted percentage of contamination [default:25] 

    --completeness         int       accepted percentage of genome completeness [default:75]

    --min-genome-length    int       minimum accepted genome length [default: 50000] 

    --local-scheduler


  
  **4.a. Metagenome assembly with read quality control analysis followed by genome binning, bin refinement and condition_free de-replication** 

     [MetagenomeAnalysis]$ metagenome.py  dRepBins  \
                                       --pre-process-reads  ``yes``   \
                                       --dRep-method ``condition_free`` \
                                       --checkM-method ``taxonomy_wf``  \
                                       --contamination ``10`` \
                                       --completeness ``80`` \
                                       --local-scheduler

    Note: Output folder in the name of ``dReplicated_bins_condition_free will`` be created at ``$PWD/metagenome_demo_analysis/dRep_bins/dReplicated_bins_condition_free /dereplicated_genomes`` which contains the de-replicated genomes

 
  
    **4.b. Metagenome assembly with read quality control analysis followed by genome binning, bin refinement and condition_based de-replication** 

    [MetagenomeAnalysis]$ metagenome.py  dRepBins  \
                                       --pre-process-reads  ``yes``   \
                                       --dRep-method ``condition_based`` \
                                       --reference-condition ``buccal_mucosa`` \
                                       --contrast-condition  ``tongue_dorsum`` \
                                       --checkM-method ``taxonomy_wf``  \
                                       --contamination ``10`` \
                                      --completeness ``80`` \
                                       --local-scheduler

    Note: 1. Output folder in the name of ``dReplicated_bins_condition_free will`` be created at ``$PWD/metagenome_demo_analysis/dRep_bins/dReplicated_bins_condition_based/`` which contains the de-replicated genomes
 



**5. annotate de-replicated genome bins**

    Annotation of condition-based or condition-free de-replicated genome bins can be done using command ``annotateGenomeBins``

    **Requirements**
    1. ``annotateGenomeBins`` must be run after condition_based or condition_free de-replication of genome bins 
    2. Pre execution of ``configureProject`` command 
    3. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
    4. If user opted for condition_based de-replication, availability of ``genomes_dRep_by_condition.lst`` file containing the list of genome bins, inside the ``$PWD/metagenome_demo_analysis/dRep_bins`` folder
    5. If user opted for condition_based de-replication, availability of ``genomes_dRep_regardless_condition.lst`` file containing the list of genome bins inside the ``$PWD/metagenome_demo_analysis/dRep_bins`` folder

    [MetagenomeAnalysis]$ metagenome.py  dRepBins <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the ilumina reads required or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

    --min-contig-length    int       minimum contig length [default 200]

    --dRep-method          str       Mandatory Parameter. choice of de-replication method. choice [condition_based, condition_free]
                    
    --local-scheduler


    Example Run
    Metagenome annotation after condition_free de-replication of genome bins** 


    [MetagenomeAnalysis]$ metagenome.py  annotateGenomeBins  \
                                       --pre-process-reads  ``yes``   \
                                       --dRep-method condition_free \                                 
                                       --local-scheduler

    Note: 1. Output folder in the name of ``prokka_condition_free`` will be created at ``$PWD/metagenome_demo_analysis/prokka_condition_free``


    Example Run
    Metagenome annotation after condition_based de-replication of genome bins** 


    [MetagenomeAnalysis]$ metagenome.py  annotateGenomeBins  \
                                       --pre-process-reads  ``yes``   \
                                       --dRep-method condition_based \                                 
                                       --local-scheduler

    Note: Output folder in the name of ``prokka_condition_based`` will be created at ``$PWD/metagenome_demo_analysis/prokka_condition_based``

|      



    6. Differential Enrichment Analysis of de-replicated genome bins

    Differential Enrichment Analysis of de-replicated genome bins can be done using command ``deaGenomeBins``

    Requirements
         1. Pre execution of ``configureProject`` command 
         2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
         3. Availability of ``metagenome_group.tsv`` file inside the ``config`` folder 

    [MetagenomeAnalysis]$ metagenome.py  deaGenomeBins <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the ilumina reads required or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

    --checkM-method        str       choice of workflow to run checkm. choice [lineage_wf, taxonomy_wf]
                                     default: taxonomy_wf

    --contamination        int       accepted percentage of contamination [default:25] 

    --completeness         int       accepted percentage of genome completeness [default:75]

    --min-genome-length    int       minimum accepted genome length [default: 50000] 

    --reference-condition  str       Mandatory Parameter: reference condition name

    --contrast-condition   str       Mandatory Parameter: contrast condition name

    --pvalcutoff           float     Default: 0.05

    --correction           str       correction method. Choose from [Bonferroni: b, 
                                                                     FDR 2-stage Benjamini-Hochberg: fdr_tsbh
                                                                     Holm: h
                                                                     Hommel: ho  
                                                                     Holm-Sidak: hs 
                                                                     Simes-Hochberg: sh
                                                                     FDR Benjamini-Yekutieli: fdr_by
                                                                     FDR 2-stage Benjamini-Krieger-Yekutieli: fdr_tsbky]
    --local-scheduler

    

B. Marker Gene based Taxonomy Profiling
----------------------------------------

    1. Taxonomy profiling from metagenomic shotgun sequencing data


    Taxonomy profiling of individual samples can be done using command ``compositionProfiling``


    **Requirements**
    1. Pre execution of ``configureProject`` command 
    2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
  
 
    [MetagenomeAnalysis]$ metagenome.py  compositionProfiling <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the raw illumina reads required or not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.


  **Composition profiling with read quality control analysis** 
  
    [MetagenomeAnalysis]$ metagenome.py  compositionProfiling  \
                                       --pre-process-reads  ``yes``                                    
                                       --local-scheduler


    Note: 1. Output folder in the name of ``profiled_samples`` will be created at ``$PWD/metagenome_demo_analysis/``, which contains the taxonomy profile of individual     samples
    

    2. Generate plots using hclust, graphlan and ampvis


    Requirements
     1. Pre execution of ``configureProject`` command 
     2. Availability of ``luigi.cfg`` file in ``parent folder``,``metagenome_group.tsv`` and ``pe_samples.lst`` files inside the ``config`` folder.

    [MetagenomeAnalysis]$ metagenome.py  compositionProfiling <arguments> --local-scheduler

    argument               type      Description

    [mandatory]
    --pre-process-reads    str       Run Quality Control Analysis of the raw illumina reads required or not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.
    
    [optional]

    --condition-column     str       Name of the condition column in metagenome_group.tsv file. [Default: conditions]
    --local-scheduler



    **Example**
    **Composition profiling with read quality control analysis** 

     [MetagenomeAnalysis]$ metagenome.py  profileTaxonomy  \
                                       --pre-process-reads  ``yes``  
                                       --condition-column ``conditions`` \                                  
                                       --local-scheduler


    Note: 1. Output folder in the name of ``figures`` will be created at ``$PWD/metagenome_demo_analysis/``, which contains the images
  
      

