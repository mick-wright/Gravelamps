'''
Lens Generation Testing

The script below uses the lens generation run output to
generate amplification factor plots for verifying the
functionality of the lens generation codes

Written by Mick Wright 2021
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"text.usetex": True, "font.family":"serif"})

#Folders containing the amplification factor data
folders = ["pointlens-wave", "pointlens-geo",
           "sislens-wave", "sislens-geo",
           "nfwlens-wave", "nfwlens-geo"]

#Read in dimensionless frequency and source position files
w = np.loadtxt("pointlens-wave/data/w.dat")
y = np.loadtxt("y.dat", dtype=str)

#Generate the legend text
legend_labels = np.core.defchararray.add("y = ", y)

#Loop through folders and generate plot for each
for folder in folders:
    freal = np.loadtxt(folder+"/data/fReal.dat")
    fimag = np.loadtxt(folder+"/data/fImag.dat")
    fabs = np.transpose(np.abs(freal + 1j*fimag))

    plt.figure(figsize=(16,10))
    plt.plot(w, fabs, linewidth=0.5)
    plt.xlabel("Dimensionless Frequency $(w)$")
    plt.ylabel("Amplification Factor $(|F|)$")
    plt.legend(legend_labels)
    plt.xscale("log")
    plt.savefig(folder+"_plot.pdf", bbox_inches="tight") 
