import os
import sys

#Pass arguments for the location of the vcf file and the destination locations
vcf_path = sys.argv[1]
dest_path  = sys.argv[2]

#File path
vcf_file = vcf_path

#Get line data from vcf file
vcf = open (vcf_file, "r")
lines = vcf.readlines()
vcf.close()


#Initialize instance variables
_samplecount = 0
sampleID = []
het_files = []
hom_files = []


for line in lines:
    #Get variant data from header
    if "#CHROM" in line:
        header = line.split("\t")
        _samplecount = len(header) - 9

        # Put locations of directories and het/hom files in lists
        for i in xrange(0, _samplecount):
            sampleID.append(header[9+i])
            os.makedirs(dest_path + sampleID[i].strip())
            het_files.append(open (dest_path + sampleID[i].strip() + "/" + sampleID[i].strip() + "hetbed.txt", "w"))
            hom_files.append(open (dest_path + sampleID[i].strip() + "/" + sampleID[i].strip() + "hombed.txt", "w"))


    #Read through each row and get field data. Is "PASS" the best field to check in each row?
    if "PASS" in line:
        edit = line.split("\t")
        chromNum = edit[0]
        Pos = int(edit[1])
        ID = edit[2]
        Ref = edit[3]
        Alt = edit[4]
        # For each variant, sort data in het/hom files in lists
        for i in xrange(0, _samplecount):
            Var = edit[9 + i]
            Separate = Var.split("|")
            if int(Separate[0]) != int(Separate[1]):
                het_files[i].write("chr" + str(chromNum) + "\t" + str(Pos-1) + "\t" + str(Pos) + "\t" + ID + "\t" + Ref + "\t" + Alt + "\n")
            else:
                hom_files[i].write("chr" + str(chromNum) + "\t" + str(Pos-1) + "\t" + str(Pos) + "\t" + ID + "\t" + Ref + "\t" + Alt + "\n")


# Close each file before finishing script
for i in xrange(0, _samplecount):
    het_files[i].close()
    hom_files[i].close()





