# script to filter indels after PASS-PASS was executed from Strelka/Mutect2 
# filter on:
#   1. No Alt evidence in the Normal sample,
#      Get allelic depth AD for normal sample. This contains 2 values, REF and ALT. ALT should be 0 in normal.
#      There should not be any indel information in the ALT allels
#   2. Strand Bias light, (pindel)
#      strand bias numbers: [REf-fwd REF-rev ALT-fwd ALT-rev]
#      nr of allelic indels on the forward and reverse on the ALT are to be at least: 2 
#      or the nr of allelic indels on the forward OR the reverse ALT is to be above 3. 
#      tumorSB[2] >=2 and tumorSB[3] >=2) or (tumorSB[2] >= 3 or tumorSB[3] >=3)
#   3. Strand Bias Heavy (real strand bias)
#      tumorSB[2] == 0 or tumorSB[3] == 0:
#      All allelic indel evidence needs to be present in both FWD and REV
#   4. Repeat filter
#      remove repeat region indels filters on RPA>9 in the ALT value INFO/RPA[1]>9 (RPA=1,2) 
#      RPA=( Nr repeats in REF, Nr repeats in ALT --> tumor sample) 
#      Using the absolute length of the insert or deletion, and filters on the maximum length of 4  and an RPA > 9
#   5. Allelic frequency filter
#      calculate Allelic frequency, tumor ALT (allelic depth:AD) / tumor Depth (DP)
#      filters on 0.3-0.7 values (should print out additional csv file with number to plot in R)
#
#   IN:
#     Python filter_shared_indels.py inputfile.PASS.vcf
#   OUT: *.NoAlt.vcf
#        *.SBtolerant.vcf
#        *.SBzero.vcf
#        *.RF.vcf
#        *.AF.vcf  

from cyvcf2 import VCF, Writer
from sys import argv

def NoALT(my_VCF,normal,outfile_NoALT):
    nn = Writer(outfile_NoALT,my_VCF,"w")
    samples=my_VCF.samples
    if samples[0] == normal:
        cl=0
    elif samples[1] == normal:
        cl=1
    else:
        print("something is wrong with the mutect output, no normal, tumor sample columns available.")

    for variants in my_VCF:
        normalAD = variants.format("AD")[cl] #allelic depth in Normal sample. [REf ALT]
        if normalAD[1] == 0: #There should be no ALT evidence in the Normal sample.
            nn.write_record(variants)
    nn.close()
    return

def SB_tolerant(inVCF,tumor,F2final):
    my_VCF = VCF(inVCF)
    nn = Writer(F2final,my_VCF,"w")
    samples=my_VCF.samples
    if samples[0] == tumor:
        cl=0
    elif samples[1] == tumor:
        cl=1
    else:
        print("something is wrong with the mutect output, no normal, tumor sample columns available.")

    for variants in my_VCF:
        tumorSB = variants.format("SB")[cl] #strand bias numbers: [REf-fwd REF-rev ALT-fwd ALT-rev]
        #print(tumorSB)
        if (tumorSB[2] >=2 and tumorSB[3] >=2) or (tumorSB[2] >= 3 or tumorSB[3] >=3):
#            print("keep: " + str(tumorSB))
            nn.write_record(variants)
    nn.close()
    return

def SB_zero(inVCF,tumor,F3final):
    my_VCF = VCF(inVCF)
    nn = Writer(F3final,my_VCF,"w")
    samples=my_VCF.samples
    if samples[0] == tumor:
        cl=0
    elif samples[1] == tumor:
        cl=1
    else:
        print("something is wrong with the mutect output, no normal, tumor sample columns available.")
    
    for variants in my_VCF:
        tumorSB = variants.format("SB")[cl] 
        if tumorSB[2] == 0 or tumorSB[3] == 0: #No supporting ALT reads on the fwd or rev strand
            print("remove: " + str(tumorSB))
            #an.write_record(variants)
        else:
            print("keep: " + str(tumorSB))      
            nn.write_record(variants)
    nn.close()
    return

def repeat(inVCF,F4final):
    my_VCF = VCF(inVCF)
    nn = Writer(F4final,my_VCF,"w")
    for variants in my_VCF:
        if variants.INFO.get("RPA"):
            ref=len(variants.REF)
            alt=len(variants.ALT[0])
            difference_length=abs(ref-alt)
            if (difference_length <= 4) and (variants.INFO.get("RPA")[1] > 9):
                print("remove: len: " + str(difference_length) + " nr repeats: " + str(variants.INFO.get("RPA")[1]))
                #failed.write_record(variants)
            else:
                nn.write_record(variants)
        else:
            nn.write_record(variants)
    nn.close()
    return

def AF(inVCF,tumor,F5final,out_csv,minvaf, maxvaf):
    my_VCF = VCF(inVCF)
    nn = Writer(F5final,my_VCF,"w")
    outfile= open(out_csv,"w")
    passed=0
    failed=0
    samples = my_VCF.samples
    if samples[0] == tumor:
        cl=0
    elif samples[1] == tumor:
        cl=1
    else:
        print("something is wrong with the mutect output, no normal, tumor sample columns available.")

    for variants in my_VCF:
        tumorAF = variants.format("AF")[cl]
        AF=tumorAF[0]
        #print(str(tumorAD_ALT) + "/" + str(tumorDP) + "=" + str(AF[0]))
#       total+=1
        if (AF > minvaf) and (AF < maxvaf):
            passed+=1
            #print(str(AF) + " passed")
            nn.write_record(variants)
        else:
            #print(str(AF) + " failed")
            failed+=1
    
        outfile.write(str(round(AF,2)))
        outfile.write("\n")
    #outfile.write(str(total) + "\n")
    outfile.write(str(passed))
    nn.close()
    return

input_vcf=argv[1]
sample_name=argv[2]
minvaf_value=argv[3]
maxvaf_value=argv[4]
minvaf = float(minvaf_value)
maxvaf = float(maxvaf_value)
F1name = sample_name + ".NoAlt.vcf"
F2name = sample_name + ".SBtolerant.vcf"
F3name = sample_name + ".SBzero.vcf"
F4name = sample_name + ".RF.vcf"
F5name = sample_name + ".AF.vcf"
F5name2 = sample_name + ".AF.csv"
inVCF = VCF(input_vcf)
header = inVCF.raw_header
n_index = header.find("##normal_sample=")
t_index = header.find("##tumor_sample=")
n_rest_header=header[n_index + len("##normal_sample="):]
t_rest_header=header[t_index + len("##tumor_sample="):]
nlindex = n_rest_header.find("\n")
tmindex = t_rest_header.find("\n")
normal = n_rest_header[:nlindex]
tumor = t_rest_header[:tmindex]
if (normal != tumor):
    NoALT(inVCF,normal,F1name)
    SB_tolerant(F1name,tumor,F2name)
    SB_zero(F2name,tumor,F3name)
    repeat(F3name,F4name)
    AF(F4name,tumor,F5name,F5name2,minvaf,maxvaf)
else:
    print("something is wrong between normal/tumor naming")
