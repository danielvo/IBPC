#! env python

import argparse
import os
import sys
import string
import subprocess

# author: Daniel Vodak (danielvo<at>ifi.uio.no)
#         Department of Tumor Biology
#         The Norwegian Radium Hospital
#         Oslo University Hospital

program_name = "IBPC value calculation"
program_version = "1.0 (2017/10/02)"
program_description = r"""A script for calculating variant support within a given sample as well as its "individual-based pool complement" ("IBPC" for short - a set of co-multiplexed samples, which do not originate from the same individual). Information necessary for each processed sample (e.g., the sample pool of origin, the associated BAM file, the variants to be investigated) need to be provided via the mandatory configuration input file. GATK tools are required for generating raw variant support information in the relevant samples. The generated information is subsequently used for calculating variant allelic fractions and IBPC support fractions."""
help_epilog = "Please consult the project GitHub page (https://github.com/danielvo/IBPC) in case of problems and/or questions."

parser = argparse.ArgumentParser(prog=program_name, description=program_description, epilog=help_epilog)
parser.add_argument("-c", "--configuration_file", type=str, required=True, help="Path to the configuration file (please consult the project GitHub page for configuration file format details: https://github.com/danielvo/IBPC)")
parser.add_argument("-g", "--GATK_jar", type=str, required=True, help="Path to GATK jar file (recommended GATK version: 3.7)")
parser.add_argument("-r", "--reference_fasta", type=str, required=True, help="Path to the reference genome fasta file (related fasta index files required by the GATK tools should also be present in the same directory)")
parser.add_argument("-j", "--java_executable", type=str, required=False, default="java", help="Path to Java 8 executable (only necessary to be specified if Java 8 is not the system's default Java version)")
parser.add_argument("-o", "--output_directory", type=str, required=True, help="Path to the intended output directory (the output directory must not exist prior to execution of this script)")
parser.add_argument("-b", "--min_base_quality", type=int, required=False, default=20, choices=range(0, 41), help="Minimum required base quality (bases with quality under the specified threshold will be disregarded during variant support calculation)")
parser.add_argument("-m", "--min_mapping_quality", type=int, required=False, default=20, choices=range(0, 41), help="Minimum required read mapping quality (bases from alignments with mapping quality under the specified threshold will be disregarded during variant support calculation)")
parser.add_argument("-p", "--parallelization_level", type=int, required=False, default=8, choices=range(1, 65), help="Number of GATK Pileup commands run in parallel")
parser.add_argument("-v", "--version", action="version", version=program_version)

args = parser.parse_args()

configuration_file = os.path.abspath(args.configuration_file)
if not (os.path.isfile(configuration_file)):
    sys.stderr.write("Specified configuration file not found. Exiting (1).\n")
    exit(1)

GATK_jar = os.path.abspath(args.GATK_jar)
if not (os.path.isfile(GATK_jar)):
    sys.stderr.write("Specified GATK jar file not found. Exiting (1).\n")
    exit(1)

reference_fasta = os.path.abspath(args.reference_fasta)
if not (os.path.isfile(reference_fasta)):
    sys.stderr.write("Specified reference fasta file not found. Exiting (1).\n")
    exit(1)

if (args.java_executable != "java"):
    java_executable = os.path.abspath(args.java_executable)
    if not (os.path.isfile(java_executable)):
        sys.stderr.write("Specified Java executable file not found. Exiting (1).\n")
        exit(1)
else:
    java_executable = args.java_executable

output_directory = os.path.abspath(args.output_directory)
if (os.path.isdir(output_directory)):
    sys.stderr.write("The specified output directory already exists. Exiting (1).\n")
    exit(1)
(output_directory_parent, output_directory_name) = os.path.split(output_directory)
if not (os.path.isdir(output_directory_parent)):
    sys.stderr.write("Parent directory of the specified output directory not found. Exiting (1).\n")
    exit(1)
else:
    os.system("mkdir " + output_directory)

pileup_commands_file_path = os.path.join(output_directory, "pileup_commands.sh")
pileup_commands_stdout_file_path = os.path.join(output_directory, "pileup_commands_stdout.log")
pileup_commands_stderr_file_path = os.path.join(output_directory, "pileup_commands_stderr.log")
IBPC_results_file_path = os.path.join(output_directory, "IBPC_values.tsv")
parallel_process_limit = args.parallelization_level


# A function for extracting variant read support information
# (alt allele occurrences and total allele occurrences
#  are counted at given variant locations).
# The read support is extracted from given pileup file ("pileup_file_path")
# at specified variant locations (listed in "variant_dictionary").
def update_read_support_values(variant_dictionary, pileup_file_path,
                               all_reads_key, alt_reads_key):
    SAMPLE_PILEUP = open(pileup_file_path, "r")

    for line in SAMPLE_PILEUP:
        line_s = line.strip().split()

        # skip lines without variant support information
        if (len(line_s) < 7):
            continue

        # identify the variant position and parse the variant pileup information
        var_pos = line_s[0] + ":" + line_s[1]
        bases = line_s[3]
        quals = line_s[4]
        read_infos = line_s[6].split(",")

        # skip the variant if it's not relevant for the given sample
        if not (var_pos in variant_dictionary):
            continue

        # look up the alt allele base
        alt = variant_dictionary[var_pos]["alt"]

        # initialize the allele occurrence counters
        alt_depth = 0
        all_depth = 0

        # track support of the alt allele as well as the total support
        # only consider alignments/bases satisfying the specified quality thresholds
        for i in range(len(bases)):
            base = bases[i]
            qual = ord(quals[i]) - 33

            read_info = read_infos[i]
            # read_offset = read_info.split("@")[1]
            MAPQ = int(read_info.split("@")[3])

            if ((MAPQ >= args.min_mapping_quality) and (qual >= args.min_base_quality)):
                if (base == alt):
                    alt_depth += 1
                all_depth += 1

        # update the allele occurrence counters
        variant_dictionary[var_pos][all_reads_key] += all_depth
        variant_dictionary[var_pos][alt_reads_key] += alt_depth

    SAMPLE_PILEUP.close()
    return


sample_info = {}
pool_info = {}
variant_pools = {}

# read in sample pool information
print "Processing information from the configuration file."
pool_config_file = open(configuration_file, "r")
line_number = 0

for line in pool_config_file:
    line_s = line.strip().split("\t")
    line_number += 1

    # ignore any header/comment lines
    if (line[0] == "#"):
        continue
    # exit immediately in case of an unexpected number of columns
    elif (len(line_s) != 8):
        sys.stderr.write("Line " + str(line_number) + " of the configuration file doesn't contain 8 tab-separated value columns. Exiting (2).\n")
        exit(2)

    pool_id = line_s[0]
    pool_size = line_s[1]
    sample_id = line_s[2]
    individual = line_s[3]
    BAM_path = line_s[4]
    VCF_path = line_s[5]
    contamination_value = line_s[6]
    contaminant_only = line_s[7]

    if (contaminant_only != "Y"):
        BAM_path = os.path.abspath(BAM_path)
        if not (os.path.isfile(BAM_path)):
            sys.stderr.write("BAM file specified on line " + str(line_number) + " of the configuration file (" + BAM_path + ") was not found. Exiting (2).\n")
            exit(2)

        VCF_path = os.path.abspath(VCF_path)
        if not (os.path.isfile(VCF_path)):
            sys.stderr.write("VCF file specified on line " + str(line_number) + " of the configuration file (" + VCF_path + ") was not found. Exiting (2).\n")
            exit(2)

    sample_info[sample_id] = {"pool_id": pool_id, "individual": individual,
                              "BAM_path": BAM_path, "VCF_path": VCF_path,
                              "contamination_value": contamination_value,
                              "contaminant_only": contaminant_only}

    if not (pool_id in pool_info):
        pool_info[pool_id] = {"pool_size": pool_size,
                              "pileup_dir": os.path.join(output_directory, pool_id + "_pileups"),
                              "variant_interval_file": os.path.join(output_directory, pool_id + "_variants.intervals"),
                              "samples": [sample_id]}
        variant_pools[pool_id] = {}

        # create a pool-specific pileup file directory
        p = subprocess.Popen(["mkdir", pool_info[pool_id]["pileup_dir"]])
        p.wait()

    else:
        pool_info[pool_id]["samples"].append(sample_id)
pool_config_file.close()


print "Step 1: merging sample variant lists into pool-specific variant interval files"
# create pool-specific variant interval files for pileup creation purposes
for sample_id in sample_info:
    if (sample_info[sample_id]["contaminant_only"] == "Y"):
        print "  Skipping contaminant-only sample: " + sample_id
    else:
        pool_id = sample_info[sample_id]["pool_id"]

        print "  Extracting variants for sample \"" + sample_id + "\" (pool \"" + pool_id + "\")"

        VCF_file = open(sample_info[sample_id]["VCF_path"], "r")

        for line in VCF_file:
            line_s = line.strip().split("\t")

            if (line[0] == "#"):
                continue

            chrom = line_s[0]
            start = line_s[1]

            variant_pools[pool_id][chrom + ":" + start + "-" + start] = ""
        VCF_file.close()

# control output with pool-wise variant counts
for pool_id in variant_pools:
    pool_variant_count = len(variant_pools[pool_id])
    pool_variant_interval_file_name = pool_info[pool_id]["variant_interval_file"]

    print "    " + str(pool_variant_count) + " unique variants in pool \"" + pool_id + "\", creating variant file \"" + pool_variant_interval_file_name + "\""

    pool_variant_interval_file = open(pool_variant_interval_file_name, "w")

    for variant in variant_pools[pool_id]:

        pool_variant_interval_file.write(variant + "\n")
    pool_variant_interval_file.close()

print "Step 2: creating pileup files"
# open the pileup-creation bash script file
pileup_commands_file = open(pileup_commands_file_path, "w")

# environment variables for the pileup-creation bash script
pileup_commands_file.write(r"""export _JAVA_OPTIONS="-XX:ParallelGCThreads=1"

GATK_JAR=""" + GATK_jar + r"""
REFERENCE_FASTA=""" + reference_fasta + r"""
JAVA_8=""" + java_executable + "\n\n")

# for each sample, create a pileup file
# based on all the variants called in the associated sample pool
sample_counter = 0
for sample_id in sample_info:

    BAM_file = sample_info[sample_id]["BAM_path"]
    pool_id = sample_info[sample_id]["pool_id"]
    interval_file = pool_info[pool_id]["variant_interval_file"]
    pileup_path = os.path.join(pool_info[pool_id]["pileup_dir"], sample_id + ".pileup")
    pileup_stdout_path = os.path.join(pool_info[pool_id]["pileup_dir"], sample_id + "_pileup_stdout.log")
    pileup_stderr_path = os.path.join(pool_info[pool_id]["pileup_dir"], sample_id + "_pileup_stderr.log")

    print "  Creating a GATK Pileup command for sample \"" + sample_id + "\" and variants from pool \"" + pool_id + "\""

    pileup_commands_file.write(r"""${JAVA_8} -jar ${GATK_JAR} \
-T Pileup \
-R ${REFERENCE_FASTA} \
-I """ + BAM_file + r""" \
-L """ + interval_file + r""" \
-verbose \
-drf DuplicateRead \
-o """ + pileup_path + r""" \
>  """ + pileup_stdout_path + r""" \
2> """ + pileup_stderr_path + r""" &""" + "\n\n")

    sample_counter += 1

    # respect the set parallelization level threshold
    if (sample_counter % parallel_process_limit == 0):

        pileup_commands_file.write(r"""wait""" + "\n\n")

# unless already present, append a 'wait' command to the end of the bash script
if (sample_counter % parallel_process_limit != 0):

    pileup_commands_file.write(r"""wait""" + "\n\n")

pileup_commands_file.close()

print "  Running the GATK pileup commands"
pileup_process = subprocess.Popen(["bash", pileup_commands_file_path, ">",
                                   pileup_commands_stdout_file_path, "2>",
                                   pileup_commands_stderr_file_path])
pileup_process.wait()

print "Step 3: calculating IBPC values"
# open the file for result output
IBPC_results_file = open(IBPC_results_file_path, "w")

# print the result file header
IBPC_results_file.write("\t".join(["SAMPLE_ID", "POOL_ID", "POOL_SIZE",
                                   "CHROMOSOME", "POSITION", "REF", "ALT",
                                   "SAMPLE_ALL_READS", "SAMPLE_VARIANT_READS",
                                   "SAMPLE_VARIANT_AF", "IBPC_ALL_READS",
                                   "IBPC_VARIANT_READS", "IBPC_VARIANT_AF",
                                   "CONPAIR_CONTAMINATION_ESTIMATE_TUMOR"]) + "\n")

# calculate all allelic fractions and IBPC support fraction values
for sample_id in sample_info:
    if (sample_info[sample_id]["contaminant_only"] == "Y"):
        print "  Skipping contaminant-only sample: " + sample_id
    else:
        pool_id = sample_info[sample_id]["pool_id"]
        pool_size = pool_info[pool_id]["pool_size"]
        VCF_path = sample_info[sample_id]["VCF_path"]
        individual = sample_info[sample_id]["individual"]
        pileup_path = os.path.join(pool_info[pool_id]["pileup_dir"], sample_id + ".pileup")
        contamination_value = sample_info[sample_id]["contamination_value"]
        pool_sample_list = pool_info[pool_id]["samples"]

        print "  Calculating IBPC values for sample \"" + sample_id + "\" (pool \"" + pool_id + "\"); other samples in the pool: " + ", ".join(pool_sample_list)

        # extract variants of the given sample (IBPC values will only be calculated for these variants)
        var_dict = {}

        VAR_LIST_FILE = open(VCF_path, "r")

        for line in VAR_LIST_FILE:
            line_s = line.strip().split("\t")

            if (line[0] == "#"):
                continue

            chr = line_s[0]
            pos = line_s[1]
            ref = line_s[3]
            alt = line_s[4]

            var_pos = chr + ":" + pos

            var_dict[var_pos] = {"chr": chr, "pos": pos,
                                 "ref": ref, "alt": alt,
                                 "sample_all_reads": 0, "sample_alt_reads": 0,
                                 "IBPC_all_reads": 0, "IBPC_alt_reads": 0}

        VAR_LIST_FILE.close()

        # skip to the next sample if no variants are to be annotated for this sample
        if (len(var_dict) == 0):
            continue

        # calculate AF values in the sample of origin
        update_read_support_values(var_dict, pileup_path, "sample_all_reads", "sample_alt_reads")

        # check for variant support in other samples in the pool,
        # unless they originate from the same individual
        for complementary_sample_id in pool_sample_list:
            complementary_individual = sample_info[complementary_sample_id]["individual"]

            if (individual == complementary_individual):
                print "    " + complementary_sample_id + " will not be used for IBPC calculation"
            else:
                print "    " + complementary_sample_id + " will be used for IBPC calculation"
                complementary_pileup_path = os.path.join(pool_info[pool_id]["pileup_dir"], complementary_sample_id + ".pileup")
                update_read_support_values(var_dict, complementary_pileup_path, "IBPC_all_reads", "IBPC_alt_reads")

        # output results for each variant of the processed sample
        for var_pos in var_dict:
            chr = var_dict[var_pos]["chr"]
            pos = var_dict[var_pos]["pos"]
            ref = var_dict[var_pos]["ref"]
            alt = var_dict[var_pos]["alt"]

            sample_all_reads = var_dict[var_pos]["sample_all_reads"]
            sample_alt_reads = var_dict[var_pos]["sample_alt_reads"]
            sample_variant_AF = 0

            if (sample_all_reads > 0):
                sample_variant_AF = round(float(sample_alt_reads)/sample_all_reads, 4)

            IBPC_all_reads = var_dict[var_pos]["IBPC_all_reads"]
            IBPC_alt_reads = var_dict[var_pos]["IBPC_alt_reads"]
            IBPC_variant_AF = 0

            if (IBPC_all_reads > 0):
                IBPC_variant_AF = round(float(IBPC_alt_reads)/IBPC_all_reads, 4)

            IBPC_results_file.write("\t".join([sample_id, pool_id, pool_size, chr, pos, ref, alt, str(sample_all_reads), str(sample_alt_reads), str(sample_variant_AF),
                                    str(IBPC_all_reads), str(IBPC_alt_reads), str(IBPC_variant_AF), contamination_value]) + "\n")

IBPC_results_file.close()
