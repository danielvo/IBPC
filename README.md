# IBPC
Calculation of IBPC read fraction variant support

The "IBPC_value_calculation_py" script requires a single mandatory argument (path to the input file) and generates multiple output files and sub-directories within the current directory. In its current form, the script requires additional parameters (mainly paths to external software files) to be set internally.

The input file can have any number of initial comment/header lines (beginning on "#"). Each data line of the input file speficies information about one of the samples that are to be processed by the script. Exactly 8 columns with the following content (and exactly in the specified order) are expected:
1) pool ID [string]: a value indicating which samples have been multiplexed together;
2) pool size [integer]: size of given pool (used in the output files, but not utilized for any calculations);
3) sample ID [string]: a value identifying a given sample;
4) individual [string]: a value indicating which samples originate from the same individual;
5) BAM file path [string]: the path to given sample's alignment file in the BAM format;
6) VCF file path [string]: the path to given sample's variant call file in the VCF format (all variants listed in the VCF file will be annotated with IBPC values);
7) contamination values [float]: given sample's contamination level (used in the output files, but not utilized for any calculations);
8) contaminant-only status ["Y"|.]: a flag indicating whether given sample should be used as a contaminantion source only (if the supplied value is "Y"), any value other than "Y" will lead to IBPC calculation for given sample's variants.
  
Created output sub-directories/files:
- a sub-directory per sample pool (i.e., per unique value in column 1 of the input file) into which pileup files of all associated samples will be placed;
- a GATK-compatible interval file per sample pool specifying all variants of a given pool that are to be included in the IBPC calculations;
- bash script file "pielup_commands.sh", which is used for running all GATK pileup commands necessary for the analysis;
- results file "IBPC_values.tsv" with IBPC-annotation information for all processed variants (includes variant-level allelic fraction and IBPC read fraction information as well as information about the associated sample and sample pool).

The code has been run with Python 2.6.6.
