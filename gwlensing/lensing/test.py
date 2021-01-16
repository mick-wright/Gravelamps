import gwlensing.lensing
import numpy as np 
import time 

w = np.loadtxt("w.dat") 
y = np.loadtxt("y.dat") 

tic = time.perf_counter()
f = gwlensing.lensing.point_lens.generate_amplification_factor_matrix(w,y)
toc = time.perf_counter() 

print(f"Took {toc - tic:0.4f} seconds") 
