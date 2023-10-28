def mu2(S_ab):
	import numpy as np
	S_abt = np.transpose(S_ab,(0,1,2,4,3))
	tr = np.trace(S_ab, axis1=3,axis2=4)**2
	mat2 = np.matmul(S_ab, S_abt)
	tr2 = np.trace(mat2, axis1=3, axis2=4)
	return (tr-tr2)*0.5

