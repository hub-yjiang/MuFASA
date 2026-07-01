from cyvcf2 import VCF, Writer
from sys import argv

input_vcf=argv[1]
sample_name=argv[2]
all_chroms=['chr' + chrom for chrom in argv[3].split(',')]
inVCF = VCF(input_vcf)
header = inVCF.raw_header
samples = inVCF.samples
outfile = sample_name + ".Manta.filtered.vcf"
outbed = open(sample_name + ".Manta.regions.bed","w") #bed file is made for validation with gridss
nn = Writer(outfile, inVCF, "w")

alt_in_normal_count = 0
# sample_name is tumor, identify the correct column in vcf samples = ['sample name',' sample name']

if sample_name in samples[0]:
    # if our tumor sample id is indeed samples[0], samples[1] will than be normal
    norm = 1
else:
    # our tumor sample id does not match the first value [0], so this is the normal sample value column in the vcf 
    norm = 0

for variant in inVCF:
    chrom = variant.CHROM
    pos = variant.POS
    if chrom in all_chroms:
        #print(chrom)
        normalPR = variant.format("PR")[norm] #paired-read information in normal sample. [Ref,ALT]
        if variant.format("SR") is not None:
            normalSR = variant.format("SR")[norm] #split-read information in normal sample. [Ref,ALT]
        else:
            normalSR = [0,0]
        if (normalPR[1] == 0) & (normalSR[1] == 0) : #There should be no ALT evidence in the Normal sample.
            # No ALT in normal sample
            nn.write_record(variant)
            print(pos)
            outbed.write(str(chrom) + "\t" + str(pos - 10000) + "\t" + str(pos + 10000) + "\n")
        else:
            alt_in_normal_count+=1
            #print(variant)
print("nr of ALT evidence in normal for:" + str(sample_name) + "  " + str(alt_in_normal_count))
nn.close()
outbed.close()
