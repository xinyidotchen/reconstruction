#Sig= 10 # smoothing scale
def G(k,smooth_scale):
	Sig=smooth_scale
	import numpy as np
	return np.exp(-0.5*k**2.*Sig**2.)
