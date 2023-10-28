def one_iteration(L,N,delta_r0,delta_L_current,delta_k,ksize,kxyz,associated_k,smooth_scale,Omega_m,b,f,weight,space,Cani):
	import gc
	def print_memory_usage():
	  import resource
	  print("Memory used:",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0)
        
	import smoothing_ani as sm
	import mu2 as mu2
	from det_I_ab import det_I_ab
	import numpy as np
	import timeit
	from scipy.interpolate import RegularGridInterpolator

	start = timeit.default_timer()
	print("Timer started....")
	
	# S_1 (k)
	print_memory_usage()

	S_l_kx=np.zeros(ksize,dtype=np.complex128)
	S_l_ky=np.zeros(ksize,dtype=np.complex128)
	S_l_kz=np.zeros(ksize,dtype=np.complex128)
#	G_element = sm.G(associated_k,smooth_scale)
	G_element=sm.G(kxyz,smooth_scale,Cani)
	S_l_kx = 1j*kxyz[0] * delta_k*G_element/associated_k**2
	S_l_ky = 1j*kxyz[1] * delta_k*G_element/associated_k**2
	S_l_kz = 1j*kxyz[2] * delta_k*G_element/associated_k**2
	S_l_kx[0,0,0] = 0.0
	S_l_ky[0,0,0] = 0.0
	S_l_kz[0,0,0] = 0.0

	print('end S_1_k')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	# S_1
	S_l_x=np.fft.irfftn(S_l_kx)/(L/N)**3.
	S_l_y=np.fft.irfftn(S_l_ky)/(L/N)**3.
	S_l_z=np.fft.irfftn(S_l_kz)/(L/N)**3.
	print("S_1x mean=",np.mean(S_l_x.flatten()))
	print("S_1x std=",np.std((S_l_x).flatten()))
	print("S_1x max=",np.max(S_l_x.flatten()))
	print("S_1x min=",np.min(S_l_x.flatten()))

	print("S_1y mean=",np.mean((S_l_y).flatten()))
	print("S_1y std=",np.std((S_l_y).flatten()))
	print("S_1y max=",np.max(S_l_y.flatten()))
	print("S_1y min=",np.min(S_l_y.flatten()))

	print("S_1z mean=",np.mean((S_l_z).flatten()))
	print("S_1z std=",np.std((S_l_z).flatten()))
	print("S_1z max=",np.max(S_l_z.flatten()))
	print("S_1z min=",np.min(S_l_z.flatten()))


	print('end S_l_1')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	# S_1ab (k)
	S_l_kxx=np.zeros(ksize,dtype=np.complex128)
	S_l_kxy=np.zeros(ksize,dtype=np.complex128)
	S_l_kxz=np.zeros(ksize,dtype=np.complex128)

	S_l_kyx=np.zeros(ksize,dtype=np.complex128)
	S_l_kyy=np.zeros(ksize,dtype=np.complex128)
	S_l_kyz=np.zeros(ksize,dtype=np.complex128)

	S_l_kzx=np.zeros(ksize,dtype=np.complex128)
	S_l_kzy=np.zeros(ksize,dtype=np.complex128)
	S_l_kzz=np.zeros(ksize,dtype=np.complex128)

	S_l_kxx=1j*kxyz[0]*S_l_kx
	S_l_kxy=1j*kxyz[1]*S_l_kx
	S_l_kxz=1j*kxyz[2]*S_l_kx

	S_l_kyx=1j*kxyz[0]*S_l_ky
	S_l_kyy=1j*kxyz[1]*S_l_ky
	S_l_kyz=1j*kxyz[2]*S_l_ky

	S_l_kzx=1j*kxyz[0]*S_l_kz
	S_l_kzy=1j*kxyz[1]*S_l_kz
	S_l_kzz=1j*kxyz[2]*S_l_kz

	print('end S_1_kab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	# S_1ab
	S_l_xx=np.fft.irfftn(S_l_kxx)/(L/N)**3.
	S_l_xy=np.fft.irfftn(S_l_kxy)/(L/N)**3.
	S_l_xz=np.fft.irfftn(S_l_kxz)/(L/N)**3.

	S_l_yx=np.fft.irfftn(S_l_kyx)/(L/N)**3.
	S_l_yy=np.fft.irfftn(S_l_kyy)/(L/N)**3.
	S_l_yz=np.fft.irfftn(S_l_kyz)/(L/N)**3.

	S_l_zx=np.fft.irfftn(S_l_kzx)/(L/N)**3.
	S_l_zy=np.fft.irfftn(S_l_kzy)/(L/N)**3.
	S_l_zz=np.fft.irfftn(S_l_kzz)/(L/N)**3.


	print("S_l_xx mean=",np.mean(S_l_xx))
	print("S_l_xx std=",np.std(S_l_xx))
	print("S_l_xx max=",np.max(S_l_xx))
	print("S_l_xx min=",np.min(S_l_xx))
	
	print("S_l_xy mean=",np.mean(S_l_xy))
	print("S_l_xy std=",np.std(S_l_xy))
	print("S_l_xy max=",np.max(S_l_xy))
	print("S_l_xy min=",np.min(S_l_xy))
	
	print("S_l_xz mean=",np.mean(S_l_xz))
	print("S_l_xz std=",np.std(S_l_xz))
	print("S_l_xz max=",np.max(S_l_xz))
	print("S_l_xz min=",np.min(S_l_xz))

	print("S_l_yx mean=",np.mean(S_l_yx))
	print("S_l_yx std=",np.std(S_l_yx))
	print("S_l_yx max=",np.max(S_l_yx))
	print("S_l_yx min=",np.min(S_l_yx))

	print("S_l_yy mean=",np.mean(S_l_yy))
	print("S_l_yy std=",np.std(S_l_yy))
	print("S_l_yy max=",np.max(S_l_yy))
	print("S_l_yy min=",np.min(S_l_yy))

	print("S_l_yz mean=",np.mean(S_l_yz))
	print("S_l_yz std=",np.std(S_l_yz))
	print("S_l_yz max=",np.max(S_l_yz))
	print("S_l_yz min=",np.min(S_l_yz))

	print("S_l_zx mean=",np.mean(S_l_zx))
	print("S_l_zx std=",np.std(S_l_zx))
	print("S_l_zx max=",np.max(S_l_zx))
	print("S_l_zx min=",np.min(S_l_zx))

	print("S_l_zy mean=",np.mean(S_l_zy))
	print("S_l_zy std=",np.std(S_l_zy))
	print("S_l_zy max=",np.max(S_l_zy))
	print("S_l_zy min=",np.min(S_l_zy))

	print("S_l_zz mean=",np.mean(S_l_zz))
	print("S_l_zz std=",np.std(S_l_zz))
	print("S_l_zz max=",np.max(S_l_zz))
	print("S_l_zz min=",np.min(S_l_zz))




	print('end S_1_ab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	S_ab=np.array([[S_l_xx,S_l_xy,S_l_xz],[S_l_yx,S_l_yy,S_l_yz],[S_l_zx,S_l_zy,S_l_zz]])
	S_ab=np.transpose(S_ab, (2,3,4,0,1))
	tr_S_l_ab = np.trace(S_ab, axis1=3,axis2=4)
	print('tr_S_l_ab mean=',np.mean(tr_S_l_ab))
	print('tr_S_l_ab std=',np.std(tr_S_l_ab))
	print('tr_S_l_ab max=',np.max(tr_S_l_ab))
	print('tr_S_l_ab min=',np.min(tr_S_l_ab))

	# Clean up memory
	del S_l_kx
	del S_l_ky 
	del S_l_kz 
	del S_l_xx
	del S_l_xy
	del S_l_xz
	del S_l_yx
	del S_l_yy
	del S_l_yz
	del S_l_zx
	del S_l_zy
	del S_l_zz
	del S_l_kxx
	del S_l_kxy
	del S_l_kxz
	del S_l_kyx
	del S_l_kyy
	del S_l_kyz
	del S_l_kzx
	del S_l_kzy
	del S_l_kzz
	


	print('end S_ab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()

	# smothed delta_r
	
	smo_delta_k=np.zeros(ksize,dtype=np.complex128)
	smo_delta_k=delta_k*G_element

	smo_delta_r=np.fft.irfftn(smo_delta_k)/(L/N)**3.
	print('smo_delta_r mean=',np.mean(smo_delta_r))
	print('smo_delta_r std=',np.std(smo_delta_r))
	print('smo_delta_r max=',np.max(smo_delta_r))
	print('smo_delta_r min=',np.min(smo_delta_r))


	#smo_delta_r=np.fft.irfftn(delta_k)
	print('end smo_delta_r')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	del smo_delta_k

	# S_2
	div_S_2=-3./14.*Omega_m**(-1./143.)*2*mu2.mu2(S_ab)
	print("div_S_2 mean=", np.mean(div_S_2))
	print("div_S_2 std=", np.std(div_S_2))
	print("div_S_2 max=", np.max(div_S_2))
	print("div_S_2 min=", np.min(div_S_2))

	print('end div_S_2')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()

	fft_div_S_2= np.fft.rfftn(div_S_2)*(L/N)**3.
	print('end fft_div_S_2')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	del div_S_2

	S_2_fftx=np.zeros(ksize,dtype=np.complex128)
	S_2_ffty=np.zeros(ksize,dtype=np.complex128)
	S_2_fftz=np.zeros(ksize,dtype=np.complex128)

	S_2_fftx=-1j*kxyz[0]*fft_div_S_2/associated_k**2.
	S_2_ffty=-1j*kxyz[1]*fft_div_S_2/associated_k**2.
	S_2_fftz=-1j*kxyz[2]*fft_div_S_2/associated_k**2.

	S_2_fftx[0,0,0]=0
	S_2_ffty[0,0,0]=0
	S_2_fftz[0,0,0]=0
	del fft_div_S_2

	print('end S_2fft')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()

	S_2x=np.fft.irfftn(S_2_fftx)/(L/N)**3.
	S_2y=np.fft.irfftn(S_2_ffty)/(L/N)**3.
	S_2z=np.fft.irfftn(S_2_fftz)/(L/N)**3.
	del S_2_fftx
	del S_2_ffty
	del S_2_fftz


	print("S_2x mean=",np.mean(S_2x))
	print("S_2x std=",np.std((S_2x).flat))
	print("S_2x max=",np.max(S_2x))
	print("S_2x min=",np.min(S_2x))

	print("S_2y mean=",np.mean(S_2y))
	print("S_2y std=",np.std((S_2y).flat))
	print("S_2y max=",np.max(S_2y))
	print("S_2y min=",np.min(S_2y))

	print("S_2z mean=",np.mean(S_2z))
	print("S_2z std=",np.std((S_2z).flat))
	print("S_2z max=",np.max(S_2z))
	print("S_2z min=",np.min(S_2z))

	print('end S_2')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	# S_s
	S_s_x=S_l_x+S_2x
	S_s_y=S_l_y+S_2y
	S_s_z=S_l_z+S_2z

	if space=='z':
		beta=f/b
		S_s_z=S_l_z+S_2z+beta*(S_l_z+2*S_2z) #rsd
	
	# turn off S_2
#	S_s_x=S_l_x
#	S_s_y=S_l_y
#	S_s_z=S_l_z+f*(S_l_z)

	print("S_s_x mean=",np.mean(S_s_x))
	print("S_s_x std=",np.std((S_s_x).flat))
	print("S_s_x max=",np.max(S_s_x))
	print("S_s_x min=",np.min(S_s_x))

	print("S_s_y mean=",np.mean(S_s_y))
	print("S_s_y std=",np.std((S_s_y).flat))
	print("S_s_y max=",np.max(S_s_y))
	print("S_s_y min=",np.min(S_s_y))

	print("S_s_z mean=",np.mean(S_s_z))
	print("S_s_z std=",np.std((S_s_z).flat))
	print("S_s_z max=",np.max(S_s_z))
	print("S_s_z min=",np.min(S_s_z))



	#Clean up
	#del S_l_x
	#del S_l_y
	#del S_l_z
	del S_2x
	del S_2y
	del S_2z

	S_s_kx= np.fft.rfftn(S_s_x)*(L/N)**3.
	S_s_ky= np.fft.rfftn(S_s_y)*(L/N)**3.
	S_s_kz= np.fft.rfftn(S_s_z)*(L/N)**3.
	#del S_s_x
	#del S_s_y
	#del S_s_z

	print('end S_s_k')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	# S_sab (k)
	S_s_kxx=np.zeros(ksize,dtype=np.complex128)
	S_s_kxy=np.zeros(ksize,dtype=np.complex128)
	S_s_kxz=np.zeros(ksize,dtype=np.complex128)

	S_s_kyx=np.zeros(ksize,dtype=np.complex128)
	S_s_kyy=np.zeros(ksize,dtype=np.complex128)
	S_s_kyz=np.zeros(ksize,dtype=np.complex128)

	S_s_kzx=np.zeros(ksize,dtype=np.complex128)
	S_s_kzy=np.zeros(ksize,dtype=np.complex128)
	S_s_kzz=np.zeros(ksize,dtype=np.complex128)

	S_s_kxx=1j*kxyz[0]*S_s_kx
	S_s_kxy=1j*kxyz[1]*S_s_kx
	S_s_kxz=1j*kxyz[2]*S_s_kx

	S_s_kyx=1j*kxyz[0]*S_s_ky
	S_s_kyy=1j*kxyz[1]*S_s_ky
	S_s_kyz=1j*kxyz[2]*S_s_ky

	S_s_kzx=1j*kxyz[0]*S_s_kz
	S_s_kzy=1j*kxyz[1]*S_s_kz
	S_s_kzz=1j*kxyz[2]*S_s_kz

	print('end S_s_kab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	# S_sab
	S_s_xx=np.fft.irfftn(S_s_kxx)/(L/N)**3.
	S_s_xy=np.fft.irfftn(S_s_kxy)/(L/N)**3.
	S_s_xz=np.fft.irfftn(S_s_kxz)/(L/N)**3.

	S_s_yx=np.fft.irfftn(S_s_kyx)/(L/N)**3.
	S_s_yy=np.fft.irfftn(S_s_kyy)/(L/N)**3.
	S_s_yz=np.fft.irfftn(S_s_kyz)/(L/N)**3.

	S_s_zx=np.fft.irfftn(S_s_kzx)/(L/N)**3.
	S_s_zy=np.fft.irfftn(S_s_kzy)/(L/N)**3.
	S_s_zz=np.fft.irfftn(S_s_kzz)/(L/N)**3.


	print("S_s_xx mean=",np.mean(S_s_xx))
	print("S_s_xx std=",np.std(S_s_xx))
	print("S_s_xx max=",np.max(S_s_xx))
	print("S_s_xx min=",np.min(S_s_xx))

	print("S_s_xy mean=",np.mean(S_s_xy))
	print("S_s_xy std=",np.std(S_s_xy))
	print("S_s_xy max=",np.max(S_s_xy))
	print("S_s_xy min=",np.min(S_s_xy))

	print("S_s_xz mean=",np.mean(S_s_xz))
	print("S_s_xz std=",np.std(S_s_xz))
	print("S_s_xz max=",np.max(S_s_xz))
	print("S_s_xz min=",np.min(S_s_xz))

	print("S_s_yx mean=",np.mean(S_s_yx))
	print("S_s_yx std=",np.std(S_s_yx))
	print("S_s_yx max=",np.max(S_s_yx))
	print("S_s_yx min=",np.min(S_s_yx))

	print("S_s_yy mean=",np.mean(S_s_yy))
	print("S_s_yy std=",np.std(S_s_yy))
	print("S_s_yy max=",np.max(S_s_yy))
	print("S_s_yy min=",np.min(S_s_yy))

	print("S_s_yz mean=",np.mean(S_s_yz))
	print("S_s_yz std=",np.std(S_s_yz))
	print("S_s_yz max=",np.max(S_s_yz))
	print("S_s_yz min=",np.min(S_s_yz))

	print("S_s_zx mean=",np.mean(S_s_zx))
	print("S_s_zx std=",np.std(S_s_zx))
	print("S_s_zx max=",np.max(S_s_zx))
	print("S_s_zx min=",np.min(S_s_zx))

	print("S_s_zy mean=",np.mean(S_s_zy))
	print("S_s_zy std=",np.std(S_s_zy))
	print("S_s_zy max=",np.max(S_s_zy))
	print("S_s_zy min=",np.min(S_s_zy))

	print("S_s_zz mean=",np.mean(S_s_zz))
	print("S_s_zz std=",np.std(S_s_zz))
	print("S_s_zz max=",np.max(S_s_zz))
	print("S_s_zz min=",np.min(S_s_zz))




	print('end S_s_ab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	del S_s_kxx
	del S_s_kxy
	del S_s_kxz
	del S_s_kyx
	del S_s_kyy
	del S_s_kyz
	del S_s_kzx
	del S_s_kzy
	del S_s_kzz
	
	S_s_Ixx=S_s_xx+1
	S_s_Iyy=S_s_yy+1
	S_s_Izz=S_s_zz+1

	print('end S_s_Iab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()
	
	S_s_Iab=np.array([[S_s_Ixx,S_s_xy,S_s_xz],[S_s_yx,S_s_Iyy,S_s_yz],[S_s_zx,S_s_zy,S_s_Izz]])
	S_s_Iab=np.transpose(S_s_Iab, (2,3,4,0,1))	
	

	S_s_ab=np.array([[S_s_xx,S_s_xy,S_s_xz],[S_s_yx,S_s_yy,S_s_yz],[S_s_zx,S_s_zy,S_s_zz]])
	S_s_ab=np.transpose(S_s_ab, (2,3,4,0,1))
#	det_S_s_ab=det_I_ab(S_s_ab)
#	mu2_S_s_ab=mu2.mu2(S_s_ab)
	
	tr_S_s_ab = np.trace(S_s_ab, axis1=3,axis2=4)
	print('tr_S_s_ab mean=',np.mean(tr_S_s_ab))
#	print('mu2(S_s_ab) mean=',np.mean(mu2_S_s_ab))
#	print('det(S_s_ab) mean=',np.mean(det_S_s_ab))
	del S_s_xx
	del S_s_yy
	del S_s_zz
	del S_s_Ixx
	del S_s_Iyy
	del S_s_Izz
	del S_s_xy
	del S_s_xz
	del S_s_yx
	del S_s_yz
	del S_s_zx
	del S_s_zy

	print('end S_s_Iab_matrix')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()	

	det_S_s_Iab=det_I_ab(S_s_Iab)
	#det_S_s_Iab=1

	print("det(S_s_Iab) mean=",np.mean(det_S_s_Iab))
	print("det(S_s_Iab) std=",np.std(det_S_s_Iab))
	print("det(S_s_Iab) max=",np.max(det_S_s_Iab))
	print("det(S_s_Iab) min=",np.min(det_S_s_Iab))

	del S_s_Iab

	print('end det_S_s_Iab')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()

	xgrid=np.linspace(0,N-1,N)
	ygrid=np.linspace(0,N-1,N)
	zgrid=np.linspace(0,N-1,N)

	inter_fun= RegularGridInterpolator((xgrid, ygrid, zgrid),delta_r0,bounds_error=False,fill_value=None)
	qmesh=np.meshgrid(xgrid,ygrid,zgrid,indexing='ij')
	qmesh=np.transpose(qmesh, (1,2,3,0))
	rx=qmesh[:,:,:,0]+S_s_x/(L/N)
	ry=qmesh[:,:,:,1]+S_s_y/(L/N)
	rz=qmesh[:,:,:,2]+S_s_z/(L/N)
	def periodic(ng, ix):
		ww = np.nonzero(ix>=(ng))
		ix[ww]-=ng
		return ix

        #periodic boundry condition for ixm=-1: -1 -> grid size-1
	def periodic2(ng,ixm):
		ww=np.nonzero((ixm)<(0))
		ixm[ww]+=(ng)
		return ixm
	rx=rx.ravel()
	rx=periodic(N,rx)
	rx=periodic2(N,rx)
	ry=ry.ravel()
	ry=periodic(N,ry)
	ry=periodic2(N,ry)
	rz=rz.ravel()
	rz=periodic(N,rz)
	rz=periodic2(N,rz)

#	print(S_s_x[0,511,10]/(L/N))
#	print(S_s_y[0,511,10]/(L/N))
#	print(S_s_z[0,511,10]/(L/N))
#	print(rx[0*N**2+511*N+10])
#	print(ry[0*N**2+511*N+10])
#	print(rz[0*N**2+511*N+10])
#	print(inter_fun([rx[0*N**2+511*N+10],ry[0*N**2+511*N+10],rz[0*N**2+511*N+10]]))
#	print(delta_r0[0,511,10])
	delta_r0_q=inter_fun(np.column_stack((rx,ry,rz)))
	
	#q_pos=np.vstack([qx[0].ravel(),qy[1].ravel(),qz[2].ravel()])
	#delta_r0_q=inter_fun(np.column_stack((q_pos[0],q_pos[1],q_pos[2])))
	delta_r0_q=np.reshape(delta_r0_q,(N,N,N))
	
	delta_L=det_S_s_Iab*(1+delta_r0_q)-1+smo_delta_r
	#delta_L=det_S_s_Iab*(1+delta_r0)-1-tr_S_l_ab
	print("delta_r0 mean=",np.mean(delta_r0))

	print("delta_r0_q mean=",np.mean(delta_r0_q))
	print("delta_r0_q std=",np.std(delta_r0_q))
	print("delta_r0_q max=",np.max(delta_r0_q))
	print("delta_r0_q min=",np.min(delta_r0_q))


	print("delta r mean before weight=",np.mean(delta_L))
	print("delta r std before weight=",np.std(delta_L))
	print("delta r max before weight=",np.max(delta_L))
	print("delta r min before weight=",np.min(delta_L))

	r_con=np.sum((delta_L-delta_L_current)**2.)/np.sum(delta_r0**2.)
	print("r_con=",r_con)

	print("weight=",weight)
	delta_L=weight*delta_L+(1-weight)*delta_L_current
	print("delta r mean before subtract=",np.mean(delta_L))
	print("delta r std before subtract=",np.std(delta_L))
	print("delta r max before subtract=",np.max(delta_L))
	print("delta r min before subtract=",np.min(delta_L))

	delta_L=delta_L-np.mean(delta_L)
	print("delta r mean=",np.mean(delta_L))
	print("delta r std=",np.std(delta_L))
	print("delta r max=",np.max(delta_L))
	print("delta r min=",np.min(delta_L))

	del det_S_s_Iab
	del smo_delta_r
	delta_k=np.array(np.fft.rfftn(delta_L))*(L/N)**3.

	print('end delta_L, delta_k')
	print("Elapsed time ", timeit.default_timer()-start)
	print_memory_usage()

	return delta_L,delta_k,r_con,S_s_x,S_s_y,S_s_z,S_l_x,S_l_y,S_l_z
