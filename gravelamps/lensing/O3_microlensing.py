from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import h5py

#Read the lookuptable
f = h5py.File('lookuptable.h5', 'r')
list(f.keys())
d11 = f['f11']
d12 = f['f12']
d21 = f['f21']
d22 = f['f22']
f.close()
D11 = []
D12 = []
D21 = []
D22 = []
for i in np.arange(201):
    D11.append(d11[i])
for i in np.arange(201):
    D12.append(d12[i])
for i in np.arange(291):
    D21.append(d21[i])
for i in np.arange(291):
    D22.append(d22[i])

D11 = np.array(D11)
D12 = np.array(D12)
D21 = np.array(D21)
D22 = np.array(D22)

d1 = D11 + 1j*D12
d2 = D21 + 1j*D22

def amplification_isolated_pm(dimensionless_frequency_array, y):        
    # Calculate Amplification functions using isolated point mass.
    mup = 0.5 + ((2+y**2) / (2*y*np.sqrt(4+y**2)))
    mum = 0.5 - ((2+y**2) / (2*y*np.sqrt(4+y**2)))
    del_T = y*np.sqrt(4+y**2)/2 + np.log((np.sqrt(4+y**2)+y)/(np.sqrt(4+y**2)-y))
    jstep = 8.0*np.pi*4.93E-6*0.1*(1+1)*1000
    w = dimensionless_frequency_array
    if max(w) <= 246.0:
        if y > 0.3:
            istep = 0.01
            ii = int((y-0.1)/istep)
            j = (w/jstep).astype(int)
            A_F_lookup = (abs(mup)**(0.5))-1j*np.exp(1j*w*del_T)*(abs(mum)**(0.5)) + (1.0/(istep*jstep)*(((ii+1)*istep+0.1-y)*((j+1)*jstep-w)*d2[ii][j]+((ii+1)*istep+0.1-y)*(w-j*jstep)*d2[ii][j+1]+(y-ii*istep-0.1)*((j+1)*jstep-w)*d2[ii+1][j]+(y-ii*istep-0.1)*(w-j*jstep)*d2[ii+1][j+1]))
        else:
            istep =0.001
            ii = int((y-0.1)/istep)
            j = (w/jstep).astype(int)
            A_F_lookup = (abs(mup)**(0.5))-1j*np.exp(1j*w*del_T)*(abs(mum)**(0.5)) +(1.0/(istep*jstep)*(((ii+1)*istep+0.1-y)*((j+1)*jstep-w)*d1[ii][j]+((ii+1)*istep+0.1-y)*(w-j*jstep)*d1[ii][j+1]+(y-ii*istep-0.1)*((j+1)*jstep-w)*d1[ii+1][j]+(y-ii*istep-0.1)*(w-j*jstep)*d1[ii+1][j+1]))
            return(A_F_lookup)
    else:
        A_F_lookup_g =[]
        A_F_lookup_w =[]
        for a in w:
            if abs(a) > 246.0:
                A_F_lookup_geo = (abs(mup)**(0.5))-1j*np.exp(1j*a*del_T)*(abs(mum)**(0.5))
                A_F_lookup_g.append(A_F_lookup_geo)
            else:
                if y > 0.3:
                    istep = 0.01
                    ii = int((y-0.1)/istep)
                    j = (a/jstep).astype(int)
                    A_F_lookup = (abs(mup)**(0.5))-1j*np.exp(1j*a*del_T)*(abs(mum)**(0.5)) + (1.0/(istep*jstep)*(((ii+1)*istep+0.1-y)*((j+1)*jstep-a)*d2[ii][j]+((ii+1)*istep+0.1-y)*(a-j*jstep)*d2[ii][j+1]+(y-ii*istep-0.1)*((j+1)*jstep-a)*d2[ii+1][j]+(y-ii*istep-0.1)*(a-j*jstep)*d2[ii+1][j+1]))
                    A_F_lookup_w.append(A_F_lookup)
                else:
                    istep =0.001
                    ii = int((y-0.1)/istep)
                    j = (a/jstep).astype(int)
                    A_F_lookup = (abs(mup)**(0.5))-1j*np.exp(1j*a*del_T)*(abs(mum)**(0.5)) +(1.0/(istep*jstep)*(((ii+1)*istep+0.1-y)*((j+1)*jstep-a)*d1[ii][j]+((ii+1)*istep+0.1-y)*(a-j*jstep)*d1[ii][j+1]+(y-ii*istep-0.1)*((j+1)*jstep-a)*d1[ii+1][j]+(y-ii*istep-0.1)*(a-j*jstep)*d1[ii+1][j+1]))
                    A_F_lookup_w.append(A_F_lookup)
        A_F_lookup_t = A_F_lookup_w + A_F_lookup_g
        A_F_lookup_t = np.array(A_F_lookup_t)
            return(A_F_lookup_t)
