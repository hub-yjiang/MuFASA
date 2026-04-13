import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from cyvcf2 import VCF
from sys import argv
from matplotlib.backends.backend_pdf import PdfPages
infile = argv[1]  # .vcf input file 
color = argv[2] # graph color
Snvfilter = argv[3] # "FiNGS/PoN/Final"
sampleID = argv[4] # sampleID
myPDF = sampleID + "_" + Snvfilter + "_VAFplot.pdf"
my_VCF = VCF(infile)
header = my_VCF.raw_header
index = header.find("##tumor_sample=")
rest_header=header[index + len("##tumor_sample="):]
nlindex = rest_header.find("\n")
tumor = rest_header[:nlindex]
samples=my_VCF.samples
VAFs = []
if samples[0] == tumor:
    cl=0
elif samples[1] == tumor:
    cl=1
else:
    print("something is wrong with the mutect output, no normal, tumor sample columns available.")

for variants in my_VCF:
    tumorAF = variants.format("AF")[cl]
    AF=tumorAF[0]
    VAFs.append(AF)
    
nrVAFs = len(VAFs)
label = "SNVs after " + Snvfilter + ": " + str(nrVAFs)
n, bins, patches = plt.hist(VAFs, bins=50, range=(0, 1), color=color, label=label)
x_ticks = np.arange(0, 1.1, 0.1)  # 0.1 step size from 0 to 1
y_ticks = np.arange(0, max(n) + 200, 100)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.legend(loc="upper right")
plt.tight_layout()
#plt.show()
with PdfPages(myPDF) as pdf:
    pdf.savefig()
