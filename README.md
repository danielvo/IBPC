# Calculation of IBPC read fraction variant support

The "IBPC_value_calculation.py" script calculates read fraction support of specified variants within a given sample as well as its "individual-based pool complement" ("IBPC" for short) - a set of co-multiplexed samples, which do not originate from the same individual. Information necessary for each processed sample (e.g., the sample pool of origin, the associated BAM file, the variants to be investigated) need to be provided via a mandatory configuration input file (as described below). GATK tools are used for generating raw variant support information in the relevant samples. The generated raw information subsequently serves as input for calculating variant allelic fractions and IBPC support fractions.

## Command-line arguments

The usage and meaning of available command-line options:

```
  -c CONFIGURATION_FILE, --configuration_file CONFIGURATION_FILE
                        Path to the configuration file
                        [mandatory]
  -g GATK_JAR, --GATK_jar GATK_JAR
                        Path to GATK jar file (recommended GATK version: 3.7)
                        [mandatory]
  -r REFERENCE_FASTA, --reference_fasta REFERENCE_FASTA
                        Path to the reference genome fasta file, which was
                        used for any prior read mapping and variant calling
                        (related fasta index files required by the GATK tools
                        should also be present in the same directory)
                        [mandatory]
  -j JAVA_EXECUTABLE, --java_executable JAVA_EXECUTABLE
                        Path to Java 8 executable (only necessary to be
                        specified if Java 8 is not the system's default Java
                        version)
                        [default value: "java"]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to the intended output directory (the output
                        directory must not exist prior to execution of this
                        script)
                        [mandatory]
  -b [0,40], --min_base_quality [0, 40]
                        Minimum required base quality (bases with quality
                        under the specified threshold will be disregarded
                        during variant support calculation)
                        [default value: 20]
  -m [0, 40], --min_mapping_quality [0, 40]
                        Minimum required read mapping quality (bases from
                        alignments with mapping quality under the specified
                        threshold will be disregarded during variant support
                        calculation)
                        [default value: 20]
  -p [1, 64], --parallelization_level [1, 64]
                        Number of GATK Pileup commands run in parallel
                        [default value: 8]
```
   
### Configuration file format

The configuration file can contain any number of comment/header lines (beginning on "#"), which will be ignored during the processing. Each data line is expected to specify information about one of the samples that are to be processed by the script. Exactly 8 tab-separated columns with the following content (and exactly in the specified order) are expected:

1) sample pool ID: a value indicating given sample's pool of origin (the value determines which samples have been multiplexed together and is used in output directory- and file-names);
2) sample pool size: size of the sample pool (passed to the final output, but not utilized for any calculations);
3) sample ID: a value identifying given sample;
4) individual ID: a value identifying the individual, from which the given sample originates (the value is used for determining individual-based sample pool complements);
5) BAM file path: the path to given sample's alignment file in the BAM format;
6) VCF file path: the path to given sample's variant call file in the VCF format (please note that all variants listed in the VCF file will be annotated with IBPC values);
7) contamination values: given sample's contamination level (passed to the final output, but not utilized for any calculations);
8) contaminant-only status: a flag indicating whether given sample should be used as a contaminantion source only (if the supplied value is "Y"), any value other than "Y" will lead to IBPC calculation for given sample's variants; values in columns 6) and 7) are currently unused for samples that serve only as contamination sources.

All values are treated as string, even though downstream processing can benefit from sample pool size and contamination values being specified in int- and float-like strings, respectively.

*Example of a configuration file:*
```
#sample_pool_ID  sample_pool_size sample_ID individual_ID BAM_path VCF_path contamination contaminant_only
ExAmp_Pool 12  sample_01 patient_01  /BAMs/sample_01.bam /VCFs/sample_01.vcf 0.7 .
ExAmp_Pool 12  sample_02 patient_01  /BAMs/sample_02.bam /VCFs/sample_02.vcf 1.2
(...)
ExAmp_Pool 12  sample_11 patient_07  /BAMs/sample_11.bam . . Y
ExAmp_Pool 12  sample_12 patient_08  /BAMs/sample_12.bam /VCFs/sample_12.vcf 0.8 .
BrAmp_Pool 12  sample_21 patient_12  /BAMs/sample_21.bam /VCFs/sample_21.vcf 0.0 .
BrAmp_Pool 12  sample_22 patient_13  /BAMs/sample_22.bam . . Y
(...)
BrAmp_Pool 12  sample_31 patient_14  /BAMs/sample_31.bam /VCFs/sample_31.vcf 0.1 .
BrAmp_Pool 12  sample_32 patient_14  /BAMs/sample_32.bam /VCFs/sample_32.vcf 0.1 .
```

## Created output sub-directories and files
The following output is created in the specified output directory:

- a subdirectory for each sample pool (i.e., for each unique value in column 1 of the input configuration file); every subdirectory should contain GATK-generated pileup results/logs for each sample from a given pool (with the exception of samples marked as contaminant-only);
- a GATK-compatible interval file per sample pool (includes interval information for variants present in any sample of the pool, with the exception of samples marked as contaminant-only; the intervals are used by GATK during pileup generation);
- bash script file called "pileup_commands.sh", which is used for running all GATK pileup commands necessary for the analysis;
- results file called "IBPC_values.tsv" with IBPC-annotation information for all processed variants (includes per-sample variant-level allelic fraction and IBPC read fraction information, as well as information about the associated sample and sample pool).

The code has been run with Python 2.6.6.

## A basic plot based on the script's output
Example of a basic plot for visualizing the IBPC calculation results:

```
library(ggplot2)
library(cowplot)
library(plyr)

vars <- read.table(<path to the IBPC_values.tsv results file>, header = T, sep="\t")

ggplot(vars, aes(SAMPLE_VARIANT_AF, IBPC_VARIANT_AF)) + 
    geom_point(col = "black", size=3, alpha=0.50) + 
    geom_point(col = "lightblue4", size=2.5, alpha=0.50) + 
    geom_point(col="black", size=0.3) + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1.19)) + 
    facet_grid(POOL_ID~.) + 
    ylab("IBPC read fraction variant support") + 
    xlab("Sample variant AF") + 
    ggtitle("Sample variant AF and IBPC read fraction variant support\ndistributions (pool-wise)") + 
    scale_y_continuous(breaks = 0.2*c(0:5)) + 
    scale_x_continuous(breaks = 0.1*c(0:10)) + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
	  background_grid(major = "xy", minor = "xy")
```

[*example output plot for results based on two 12-sample pools*](./Basic_plot_example.png)
