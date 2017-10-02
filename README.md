# Calculation of IBPC read fraction variant support

The "IBPC_value_calculation.py" script calculates read fraction support for specified variants within a given sample as well as its "individual-based pool complement" ("IBPC" for short) - a set of co-multiplexed samples, which do not originate from the same individual. Information necessary for each processed sample (e.g., the sample pool of origin, the associated BAM file, the variants to be investigated) need to be provided via a mandatory configuration input file (as described below). GATK tools are used for generating raw variant support information in the relevant samples. The generated raw information subsequently serves as input for calculating variant allelic fractions and IBPC support fractions.

## Command-line arguments

The usage and meaning of available command-line options.

```
  -c CONFIGURATION_FILE, --configuration_file CONFIGURATION_FILE
                        Path to the configuration file
                        [mandatory]
  -g GATK_JAR, --GATK_jar GATK_JAR
                        Path to GATK jar file (recommended GATK version: 3.7)
                        [mandatory]
  -r REFERENCE_FASTA, --reference_fasta REFERENCE_FASTA
                        Path to the reference genome fasta file (related fasta
                        index files required by the GATK tools should also be
                        present in the same directory)
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

The configuration file can contain any number of comment/header lines (beginning on "#"), which will be ignored during the processing. Each data line is expected to specify information about one of the samples that are to be processed by the script. Exactly 8 columns with the following content (and exactly in the specified order) are expected:

1) sample pool ID: a value indicating given sample's pool of origin (the value determines which samples have been multiplexed together and is used in output directory- and file-names);
2) sample pool size: size of the sample pool (passed to the final output, but not utilized for any calculations);
3) sample ID: a value identifying given sample;
4) individual ID: a value identifying the individual, from which the given sample originates (the value is used for determining individual-based sample pool complements);
5) BAM file path: the path to given sample's alignment file in the BAM format;
6) VCF file path: the path to given sample's variant call file in the VCF format (please note that all variants listed in the VCF file will be annotated with IBPC values);
7) contamination values: given sample's contamination level (passed to the final output, but not utilized for any calculations);
8) contaminant-only status: a flag indicating whether given sample should be used as a contaminantion source only (if the supplied value is "Y"), any value other than "Y" will lead to IBPC calculation for given sample's variants.

All values are treated as string, even though downstream processing can benefit from sample pool size and contamination values being specified in int- and float-like strings, respectively.

## Created output sub-directories and files
The following output is created in the specified outut directory:

- a sub-directory per sample pool (i.e., per unique value in column 1 of the input file), into which pileup files of all associated samples will be placed;
- a GATK-compatible interval file per sample pool, specifying all variants of a given pool that are to be included in the IBPC calculations;
- bash script file "pielup_commands.sh", which is used for running all GATK pileup commands necessary for the analysis;
- results file "IBPC_values.tsv" with IBPC-annotation information for all processed variants (includes variant-level allelic fraction and IBPC read fraction information, as well as information about the associated sample and sample pool).

The code has been run with Python 2.6.6.
