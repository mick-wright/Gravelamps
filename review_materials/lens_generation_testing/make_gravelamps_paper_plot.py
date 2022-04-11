'''
The script below generates the amplification factor plot used in the Gravelamps paper.

Written by Mick Wright 2022
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
matplotlib.rcParams.update({"text.usetex": True, "font.family":"serif", "axes.titlesize":"30"})

#Read in the dimensionless frequency and source position files
w = np.loadtxt("pointlens-wave/data/w.dat")
y = np.loadtxt("y.dat", dtype=str)

#Generate the legend text
legend_labels = np.core.defchararray.add("y = ", y)
legend_labels = np.core.defchararray.add(legend_labels, " (Wave Optics)")

legend_labels_geo = np.core.defchararray.add("y = ", y)
legend_labels_geo = np.core.defchararray.add(legend_labels_geo, " (Geometric Optics)")

legend_labels = np.concatenate((legend_labels, legend_labels_geo))

#Read in the amplification factor data
point_wave_freal = np.loadtxt("pointlens-wave/data/fReal.dat")
point_wave_fimag = np.loadtxt("pointlens-wave/data/fImag.dat")

point_geo_freal = np.loadtxt("pointlens-geo/data/fReal.dat")
point_geo_fimag = np.loadtxt("pointlens-geo/data/fImag.dat")

sis_wave_freal = np.loadtxt("sislens-wave/data/fReal.dat")
sis_wave_fimag = np.loadtxt("sislens-wave/data/fImag.dat")

sis_geo_freal = np.loadtxt("sislens-geo/data/fReal.dat")
sis_geo_fimag = np.loadtxt("sislens-geo/data/fImag.dat")

nfw_wave_freal = np.loadtxt("nfwlens-wave/data/fReal.dat")
nfw_wave_fimag = np.loadtxt("nfwlens-wave/data/fImag.dat")

nfw_geo_freal = np.loadtxt("nfwlens-geo/data/fReal.dat")
nfw_geo_fimag = np.loadtxt("nfwlens-geo/data/fImag.dat")

#Construct the full complex amplification factor
point_wave = np.transpose(np.abs(point_wave_freal + 1j*point_wave_fimag))
point_geo = np.transpose(np.abs(point_geo_freal + 1j*point_geo_fimag))
sis_wave = np.transpose(np.abs(sis_wave_freal + 1j*sis_wave_fimag))
sis_geo = np.transpose(np.abs(sis_geo_freal + 1j*sis_geo_fimag))
nfw_wave = np.transpose(np.abs(nfw_wave_freal + 1j*nfw_wave_fimag))
nfw_geo = np.transpose(np.abs(nfw_geo_freal + 1j*nfw_geo_fimag))

#Construct the full plot
fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True)
fig.set_size_inches(16,10)
fig.supxlabel("Dimenesionless Frequency $(w)$", fontsize=20)
fig.supylabel("Amplification Factor $(|F|)$", fontsize=20)

#Point Lens Plot
ax1.plot(w, point_wave, linewidth=0.5)
ax1.plot(w, point_geo, linewidth=0.5, linestyle="--")

#SIS Plot
ax2.plot(w, sis_wave, linewidth=0.5)
ax2.plot(w, sis_geo, linewidth=0.5, linestyle="--")

#NFW Plot
ax3.plot(w, nfw_wave, linewidth=0.5)
ax3.plot(w, nfw_geo, linewidth=0.5, linestyle="--")

#Axial Control
for ax in (ax1, ax2, ax3):
    ax.set(xscale="log")

for ax in fig.get_axes():
    ax.label_outer()
    handles = ax.get_lines()

fig.legend(handles, legend_labels, loc=(0.05,0.68), fontsize=20)
fig.tight_layout()

#plt.show()
plt.savefig("paper_plot.eps", bbox_inches="tight")
