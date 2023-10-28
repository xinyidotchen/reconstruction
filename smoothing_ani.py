def G(kxyz,Sig_perp,Cani):
	import numpy as np
	associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
	k_par=kxyz[2] #z direction
	k_perp=np.sqrt(kxyz[0]**2+kxyz[1]**2)
	return np.exp(-0.5*(k_perp**2+k_par**2*Cani**2)*Sig_perp**2)
