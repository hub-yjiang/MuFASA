 # this script uses Mutect2 output to filter on Allelic frequency, and user set (min/max) values.
 # input:
 #    args1: input file
 #    args2: output .tab/.csv file with ALL vaf numbers. The last number is the total number of variants that would PASS the set thresholds.
 #    args3: output .vcf with all PASSED left-over variants
 #    args4: output .tab/.csv file with PASS vaf numbers (can be used for histogram plotting later)
 #    args5: min vaf set ( 0.3 )
 #    args6: max vaf set ( 0.7 )
from cyvcf2 import VCF, Writer
from sys import argv
allfile= open(argv[2],"w")
outvcf=argv[3] 
infile=argv[1] #input file.
passfile=open(argv[4],"w")
min_vaf_value=float(argv[5])
max_vaf_value=float(argv[6])
my_VCF = VCF(infile)
out = Writer(outvcf,my_VCF,"w")
passed=0
failed=0
header = my_VCF.raw_header
index = header.find("##tumor_sample=")
rest_header=header[index + len("##tumor_sample="):]
nlindex = rest_header.find("\n")
tumor = rest_header[:nlindex]
samples=my_VCF.samples
if samples[0] == tumor:
    cl=0
    print("column 0 is tumor")
elif samples[1] == tumor:
    print("column 1 is tumor")
    cl=1
else:
    print("something is wrong with the mutect output, no normal, tumor sample columns available.")

for variants in my_VCF:
    tumorAF = variants.format("AF")[cl]
    AF=tumorAF[0]
    if (AF > min_vaf_value) and (AF < max_vaf_value):
        passed+=1
        #print(str(AF) + " passed")
        out.write_record(variants)
        passfile.write(str(round(AF,2)))
        passfile.write("\n")
    
    allfile.write(str(round(AF,2)))
    allfile.write("\n")
allfile.write(str(passed))

out.close()
my_VCF.close()
allfile.close()
passfile.close()