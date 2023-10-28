#def det_I_ab(S_Iab,n):
#	import numpy as np
#	det_I_ab=np.zeros((n,n,n),dtype='f8')
#	for ix in range(n):
#		for iy in range(n):
#			for iz in range(n):
#				mat = S_Iab[:,:,ix,iy,iz]
#				det_I_ab[ix,iy,iz] = np.linalg.det(mat)
#	return det_I_ab


def det_I_ab(S_Iab):
	import numpy as np
	det_I_ab = np.linalg.det(S_Iab)
	return det_I_ab
