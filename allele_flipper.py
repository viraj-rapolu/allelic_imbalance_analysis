"""
Argument 1 should be the location of the indexed .bam file
Argument 2 should be the location of the .bed file
Argument 3 should be a directory where the fastq file outputs will be
"""
import os
import pysam
import string
import sys

# For the command line
bam_path = sys.argv[1]
bed_path = sys.argv[2]
dest_path = sys.argv[3]

# Get .bed data
bed = open(bed_path, "r")
lines = bed.readlines()
bed.close()

# Arrays for .bed data
pos = []
chrom = []
ref = []
alt = []
i = 1
beg_pos_dict = {}
end_pos_dict = {}

# Initialize some array data
beg_pos_dict['chr1'] = 0
chrom.insert(0, 'chr1')

# Store .bed data, note indexes in dicts
for line in lines:
    edit = line.split("\t")
    pos.append(int(edit[1]))
    chrom.append(edit[0])
    ref.append(edit[4].strip())
    alt.append(edit[5].strip())
    if chrom[i] != chrom[i - 1]:
        beg_pos_dict[chrom[i]] = i
        end_pos_dict[chrom[i - 1]] = i - 1
    i = i + 1

# Create fastq files
fastq1 = open(dest_path + "fastq1.txt", "w")
fastq2 = open(dest_path + "fastq2.txt", "w")

# Create translation table
trans = string.maketrans('ATGC', 'TACG')

i = 0 # TSF ONLY START AT 0 THE FIRST ROUND
next_start_i = 0
# Sort through .bam data
bamfile = pysam.AlignmentFile(bam_path, "rb")
for read in bamfile:
    chromosome = read.reference_name
    read_start_pos = int(read.reference_start)
    read_length = int(read.infer_read_length()) - 1
    read_end_pos = int(read_start_pos + read_length)
    i = beg_pos_dict[chromosome]
    i = next_start_i
    read_list = list(read.query_sequence)
    updated_bases = "".join(read_list)

    # START TSF MODIFIED
    if chrom[i].strip() != chromosome.strip():
        i = beg_pos_dict[chromosome]
        next_start_i = i
    while (read_start_pos > pos[i]) and chromosome == chrom[i]:
        i = i + 1
    next_start_i = i - 1
    # END TSF MODIFIED

    while read_end_pos >= pos[i] and i <= end_pos_dict[chromosome]:
        position = pos[i]
        if chrom[i].strip() == chromosome.strip():
            if int(read_start_pos) <= int(position) <= int(read_end_pos):
                read_index = position - read_start_pos
                old_base = read_list[read_index]
                if ref[i] == old_base: # Cases for swapping alleles
                    read_list[read_index] = alt[i]
                elif alt[i] == old_base:
                    read_list[read_index] = ref[i]
                elif ref[i].translate(trans) == old_base:
                    read_list[read_index] == alt[i].translate(trans)
                elif alt[i].translate(trans) == old_base:
                    read_list[read_index] = ref[i].translate(trans)
                # There is another case where the base doesn't match
                updated_bases = "".join(read_list)
        i = i + 1
    # Write all reads into fastq
    if read.flag > 128:
        fastq2.write(">%s\n%s\n+\n%s\n" % (
            read.query_name, updated_bases, "".join((chr(x + 33)) for x in read.query_qualities)))
    elif read.flag < 128:
        fastq1.write(">%s\n%s\n+\n%s\n" % (
            read.query_name, updated_bases, "".join((chr(x + 33)) for x in read.query_qualities)))
    else:
        print("Program error: read flag not in range")
        sys.exit()
bamfile.close()
