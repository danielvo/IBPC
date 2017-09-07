import sys, string, subprocess
# version 0.9 (2017/09/07) 


# A function for extracting variant read support information (all reads and alt reads are counted at given locations).
# The read support is extracted from given pileup file ("pileup_file_path")
# at specified variant locations (listed in "variant_dictionary").
def update_read_support_values(variant_dictionary, pileup_file_path, all_reads_key, alt_reads_key):
  SAMPLE_PILEUP = open(pileup_file_path, "r")

  for line in SAMPLE_PILEUP:
    line_s = line.strip().split()

    # skip lines without variant support information
    if (len(line_s) < 7):
      continue

    # identify variant
    var_pos = line_s[0] + ":" + line_s[1]
    bases = line_s[3]
    quals = line_s[4]
    read_infos = line_s[6].split(",")

    # skip the variant if it's not relevant for given sample
    if not (var_pos in variant_dictionary):
      continue

    alt = variant_dictionary[var_pos]["alt"]
  
    alt_depth = 0
    all_depth = 0

    # track support of the alt allele as well as the total support
    # only consider reads satisfying set quality limits
    for i in range(len(bases)):
      base = bases[i]
      qual = ord(quals[i]) - 33

      read_info = read_infos[i]
      read_offset = read_info.split("@")[1]
      MAPQ = int(read_info.split("@")[3])

      if ((MAPQ >= 20) and (qual >= 20)):
        if (base == alt):
          alt_depth += 1
        all_depth += 1

    # update the read support counters of given variant
    variant_dictionary[var_pos][all_reads_key] += all_depth
    variant_dictionary[var_pos][alt_reads_key] += alt_depth

  SAMPLE_PILEUP.close()
  return

# input files
pool_config_file = open(sys.argv[1], "r")

# output files
pileup_command_file = open("pileup_commands.sh", "w")
IBPC_results_file   = open("IBPC_values.tsv", "w")

parallel_process_limit = 32

# environment variables for a pileup creation bash script
pileup_command_file.write(r"""export _JAVA_OPTIONS="-XX:ParallelGCThreads=1"

GATK_KEY="<path_to_GATK_licence_key>"
P_GATK="<path_to_GATK_jar_file>"
P_GENREF="<path_to_reference_genome_fasta_file>"
JAVA_8="<path_to_java_8_executable>"
""" + "\n\n")

sample_info = {}
pool_info = {}
variant_pools = {}


# read in sample pool information
for line in pool_config_file:
  line_s = line.strip().split("\t")

  if (line[0] == "#"):
    continue

  pool_id    = line_s[0]
  pool_size  = line_s[1]
  sample_id  = line_s[2]
  individual = line_s[3]
  BAM_path   = line_s[4]
  VCF_path   = line_s[5]
  contamination_value = line_s[6]
  contaminant_only    = line_s[7]

  sample_info[sample_id] = {"pool_id":pool_id, "individual":individual, "BAM_path":BAM_path, "VCF_path":VCF_path, "contamination_value":contamination_value, "contaminant_only":contaminant_only}

  if not (pool_id in pool_info):
    pool_info[pool_id] = {"pool_size":pool_size, "pileup_dir":pool_id + "_pileups", "variant_file":pool_id + "_variants.intervals", "IBPC_file":pool_id + "_IBPC_values.tsv", "samples":[sample_id]}
    variant_pools[pool_id] = {}

    # create pool-specific pileup file directories
    p = subprocess.Popen(["mkdir", pool_id + "_pileups"])
    p.wait()

  else:
    pool_info[pool_id]["samples"].append(sample_id)
pool_config_file.close()


print "Step 1: merging sample variant lists into pool-specific variant interval files"
# create pool-specific variant interval files for pileup creation purposes
for sample in sample_info:
  if (sample_info[sample]["contaminant_only"] == "Y"):
    print "  Skipping contaminant-only sample: " + sample
  else:
    pool = sample_info[sample]["pool_id"]

    print "  Extracting variants for sample \"" + sample + "\" (pool \"" + pool + "\")"

    VCF_file = open(sample_info[sample]["VCF_path"], "r")

    for line in VCF_file:
      line_s = line.strip().split("\t")

      if (line[0] == "#"):
        continue

      chrom = line_s[0]
      start = line_s[1]
  
      variant_pools[pool][chrom + ":" + start + "-" + start] = ""
    VCF_file.close()

# control output with pool-wise variant counts
for pool in variant_pools:
  pool_variant_count = len(variant_pools[pool])
  pool_variant_file_name = pool_info[pool]["variant_file"]

  print "    " + str(pool_variant_count) + " unique variants in pool \"" + pool + "\", creating variant file \"" + pool_variant_file_name + "\""

  pool_variant_file = open(pool_variant_file_name, "w")

  for variant in variant_pools[pool]:
  
    pool_variant_file.write(variant + "\n")
  pool_variant_file.close()
  

print "Step 2: creating pileup files"
# create sample-wise pileup files
sample_counter = 0
for sample in sample_info:

  BAM_file = sample_info[sample]["BAM_path"]
  pool = sample_info[sample]["pool_id"]
  interval_file = pool_info[pool]["variant_file"]
  pileup_path = pool_info[pool]["pileup_dir"] + "/" + sample + ".pileup"
  pileup_stdout_path = pool_info[pool]["pileup_dir"] + "/" + sample + "_pileup_stdout.log"
  pileup_stderr_path = pool_info[pool]["pileup_dir"] + "/" + sample + "_pileup_stderr.log"

  print "  Compiling a pileup file for pool variant locations of sample \"" + sample + "\" (pool \"" + pool + "\")"

  pileup_command_file.write(r"""${JAVA_8} -jar ${P_GATK} \
-T Pileup \
-R ${P_GENREF} \
-I """ + BAM_file + r""" \
-L """ + interval_file + r""" \
-verbose \
-drf DuplicateRead \
-o """ + pileup_path + r""" \
>  """ + pileup_stdout_path + r""" \
2> """ + pileup_stderr_path + r""" &""" + "\n\n")

  sample_counter += 1

  if (sample_counter%parallel_process_limit == 0):

    pileup_command_file.write(r"""wait""" + "\n\n")

if (sample_counter%parallel_process_limit != 0):

  pileup_command_file.write(r"""wait""" + "\n\n")

pileup_command_file.close()

p = subprocess.Popen(["bash", "pileup_commands.sh", ">", "pileup_commands_stdout.log", "2>", "pileup_commands_stderr.log"])
p.wait()


# print result file header
IBPC_results_file.write("\t".join(["SAMPLE_ID", "POOL_ID", "POOL_SIZE", "CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE_ALL_READS", "SAMPLE_VARIANT_READS", "SAMPLE_VARIANT_AF", "IBPC_ALL_READS", "IBPC_VARIANT_READS", "IBPC_VARIANT_AF", "CONPAIR_CONTAMINATION_ESTIMATE_TUMOR"]) + "\n")

print "Step 3: calculating IBPC values"
# calculate IBPC results
for sample in sample_info:
  if (sample_info[sample]["contaminant_only"] == "Y"):
    print "  Skipping contaminant-only sample: " + sample
  else:
    pool = sample_info[sample]["pool_id"]
    pool_size = pool_info[pool]["pool_size"]
    VCF_path = sample_info[sample]["VCF_path"]
    pileup_path = pool_info[pool]["pileup_dir"] + "/" + sample + ".pileup"
    contamination_value = sample_info[sample]["contamination_value"]
    pool_sample_list = pool_info[pool]["samples"]

    print "  Calculating IBPC values for sample \"" + sample + "\" (pool \"" + pool + "\") other samples in the pool: " + str(pool_sample_list)

    # extract variants of given sample (IBPC values will only be calculated for these variants)
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

      var_dict[var_pos] = {"chr":chr, "pos":pos, "ref":ref, "alt":alt, "sample_all_reads":0, "sample_alt_reads":0, "IBPC_all_reads":0, "IBPC_alt_reads":0}

    VAR_LIST_FILE.close()

    # skip to the next sample if no variants are to be annotated for this sample
    if (len(var_dict) == 0):
      continue

    # calculate AF values in the ample of origin
    update_read_support_values(var_dict, pileup_path, "sample_all_reads", "sample_alt_reads")

    individual = sample_info[sample]["individual"]

    # check for variant support in other samples in the pool, unless they originate from the same individual
    for complementary_sample in pool_sample_list:
      complementary_individual = sample_info[complementary_sample]["individual"]

      if (individual == complementary_individual):
        print "    " + complementary_sample + " will not be used for IBPC calculation"
      else:
        print "    " + complementary_sample + " will be used for IBPC calculation"
        complementary_pileup_path = pool_info[pool]["pileup_dir"] + "/" + complementary_sample + ".pileup"
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

      IBPC_results_file.write("\t".join([sample, pool, pool_size, chr, pos, ref, alt, str(sample_all_reads), str(sample_alt_reads), str(sample_variant_AF),
                                str(IBPC_all_reads), str(IBPC_alt_reads), str(IBPC_variant_AF), contamination_value]) + "\n")

IBPC_results_file.close()
