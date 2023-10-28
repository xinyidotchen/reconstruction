def reconstruction(delta_r_bf_rec,delta_k_bf_rec,delta_k_ini,mono_pk_ini,N,L,run_n,redshift,recon_params,weight,space,dirname,header2_hada,Cani,pkxibins_spec):
	import hada_one_iteration_general as one
	import power_spectrum_mono_general as ps
#	import power_spectrum_multi_general as psm
	import power_spectrum_multi_general_pkxi as psm
	import gk_rk_general as cos
	import timeit
	import numpy as np
	import smoothing as sm
	import os
#	import dis_cor_rsd as dis
	import correlation_angle as cor_angle
#	import cross_cor_angle as cross_angle
	import rk_avg_later_general as rk_avg_later
	import amplitude_recovery_general as amp



	max_iter=30
	smooth_scale=recon_params.smoothing
	Omega_m=recon_params.OmegaM(redshift)
	b=recon_params.bias
	f=recon_params.f
	
	delta_L_current=delta_r_bf_rec/b
	delta_k_current=delta_k_bf_rec/b
	delta_r0=delta_L_current

	associated_kx=np.fft.fftfreq(N,d=1)*2*np.pi*N/L # kx vector. Python gives (0,1,...)/N. k=2*pi*q/L, so multiply it by 2*pi*N/L
	associated_ky=np.fft.fftfreq(N,d=1)*2*np.pi*N/L
	associated_kz=np.fft.fftfreq(N,d=1)[0:int(N/2+1)]*2*np.pi*N/L

	ksize = (len(associated_kx), len(associated_ky), len(associated_kz))
	# Define k grid and associated k
	kxyz = np.meshgrid(associated_kx, associated_ky, associated_kz, indexing='ij')
	associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
	associated_k[0,0,0] = 1.0e-10 # changed from 1.0 to 1.0e-10 to be consistent with other associated_k[0,0,0]


	
	rcon=1
	
#	S_tru_x=np.load('displacement_tru_x_on_grid_rsd_run%i_z%i.dat.npy'%(run_n,redshift))
#	S_tru_y=np.load('displacement_tru_y_on_grid_rsd_run%i_z%i.dat.npy'%(run_n,redshift))
#	S_tru_z=np.load('displacement_tru_z_on_grid_rsd_run%i_z%i.dat.npy'%(run_n,redshift))
	for i in range(max_iter):
	
		if rcon>=0.01:
			sm=max(20/(1.2**i), smooth_scale)	
			print('start iter '+str(i))
			start=timeit.default_timer() # time it
			print('start one_itertaion')
			delta_L_current,delta_k_current,rcon,S_s_x,S_s_y,S_s_z,S_l_x,S_l_y,S_l_z=one.one_iteration(L,N,delta_r0,delta_L_current,delta_k_current,ksize,kxyz,associated_k,sm,Omega_m,b,f,weight,space,Cani)
#			np.save('delta_k_current_%i_run%i_sm%i_%i.dat'%(n,run_n,smooth_scale,i),delta_k_current)
#			np.save('delta_r_current_%i_run%i_sm%i_%i_z%i.dat'%(n,run_n,smooth_scale,i,redshift),delta_L_current)
			
			######## saving S_s for test of S_s cross-correlation of snapshots and halo ##########
#			filename='/home/xc298/scratch60/reconstruction/output/'\
#				+dirname+'/large'+'/S_s_x_space%s_ani%.1f_iter%i.dat'%(space,Cani,i)
#			os.makedirs(os.path.dirname(filename),exist_ok=True)
#			np.save(filename,S_s_x)
#			filename='/home/xc298/scratch60/reconstruction/output/'\
#				+dirname+'/large'+'/S_s_y_space%s_ani%.1f_iter%i.dat'%(space,Cani,i)
#			np.save(filename,S_s_y)
#			filename='/home/xc298/scratch60/reconstruction/output/'\
#				+dirname+'/large'+'/S_s_z_space%s_ani%.1f_iter%i.dat'%(space,Cani,i)
#			np.save(filename,S_s_z)
			###################################################################################### 

			delta_k_current=np.array(np.fft.fftn(delta_L_current))*(L/N)**3. # changed from rfftn to fftn to run analyses
			print("r_con=", rcon)
			print("mean delta_r =",np.mean(delta_L_current))
			print('start power_spectrum')
			#k_avg,r_avg,mono_pk_current,mono_cf_current,k_counts,r_counts=ps.power_spectrum(N,L,delta_k_current,pkxibins_spec)
			k_pk,r_cf=ps.power_spectrum(N,L,delta_k_current,pkxibins_spec)
			#k_pk=np.array([k_avg,mono_pk_current,k_counts]).T
			#r_cf=np.array([r_avg,mono_cf_current,r_counts]).T
			mono_pk_current=k_pk.T[2] # w/o window
			
			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/mono_pk_hada_space%s_ani%.1f_iter%i.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, k_pk, fmt=['%.12f','%.12f','%.12f','%.12f'],\
				header="Monopole power spectrum "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s\t%12s'%('k mean','P(k) w/ window','P(k) w/o window','k counts'),delimiter='\t')

			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/mono_cf_hada_space%s_ani%.1f_iter%i.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, r_cf, fmt=['%.12f','%.12f','%.12f'],\
				header="Monopole correlation function "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s'%('r mean','xi(r)','r counts'),delimiter='\t')


			k_pk0_pk2_numPl3,r_xi0_xi2=psm.power_spectrum(N,L,delta_L_current,'hada',numPl=3,pkxibins_spec=pkxibins_spec)
			k_pk0_pk2_numPl5,r_xi0_xi2=psm.power_spectrum(N,L,delta_L_current,'hada',numPl=5,pkxibins_spec=pkxibins_spec)

			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/multi_pk_hada_space%s_ani%.1f_iter%i_numPl3.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, k_pk0_pk2_numPl3, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f','%.12f'],\
				header="Mono, quad power spectrum "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s\t%12s\t%12s\t%12s'%('k mean','P0(k) w/ window','P2(k) w/ window','P0(k) w/o window','P2(k) w/o window','k counts'),delimiter='\t')

			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/multi_pk_hada_space%s_ani%.1f_iter%i_numPl5.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, k_pk0_pk2_numPl5, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f','%.12f'],\
				header="Mono, quad power spectrum "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s\t%12s\t%12s\t%12s'%('k mean','P0(k) w/ window','P2(k) w/window','P0(k) w/o window','P2(k) w/o window','k counts'),delimiter='\t')


			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/multi_cf_hada_space%s_ani%.1f_iter%i.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, r_xi0_xi2, fmt=['%.12f','%.12f','%.12f','%.12f'],\
				header="Mono, quad correlation function "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s\t%12s'%('r mean','xi0(r)','xi2(r)','r counts'),delimiter='\t')


			
			print('start cross_cor_coeff')
			k_avg,gk_current,rk_current,numer_current=cos.cross_cor_coeff(delta_k_ini,mono_pk_ini,delta_k_current,mono_pk_current,L,N,pkxibins_spec)
			k_avg,rk_avg_later_current=rk_avg_later.rk_avg_later_cross_cor_coeff(delta_k_ini,delta_k_current,L,N,pkxibins_spec)
			k_gk_rk=np.array([k_avg,gk_current,rk_current,rk_avg_later_current,numer_current]).T

			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/gk_rk_hada_space%s_ani%.1f_iter%i.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, k_gk_rk, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f'],\
				header="G(k) and r(k) "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s\t%12s\t%12s'%('k','g(k)','r(k)','r(k) avg later','numer'),delimiter='\t')

			k_avg,cor_avg1,cor_avg2,cor_avg3,pk_avg1,pk_avg2,pk_avg3=cor_angle.correlation_angle(delta_k_current,delta_k_ini,L,N,pkxibins_spec)
			cor_w_angle=np.array([k_avg,cor_avg1,cor_avg2,cor_avg3,pk_avg1,pk_avg2,pk_avg3]).T
			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/cor_angle_numer_hada_space%s_ani%.1f_iter%i.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			np.savetxt(filename, cor_w_angle, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f','%.12f','%.12f'],\
				header="cross correlation, angle (numerator) "+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s'%('k','mu<1/3 numer','1/3<mu<2/3 numer','mu>2/3 numer','mu<1/3 pk','1/2<mu<2/3 pk','mu>2/3 pk'),delimiter='\t')

			k_avg,amplitude_ratio=amp.amplitude_recovery(delta_k_ini,delta_k_current,L,N,pkxibins_spec)
			k_amp=np.array([k_avg,amplitude_ratio]).T
			filename='/home/xc298/project/reconstruction_project/output/'\
				+dirname+'/amplitude_ratio_space%s_ani%.1f_iter%i.txt'%(space,Cani,i)
			os.makedirs(os.path.dirname(filename),exist_ok=True)
			np.savetxt(filename,k_amp,fmt=['%.12f','%.12f'],\
				header='amplitude ratio'+header2_hada+"iter %i\n"%(i)+"Current smoothing: %f Mpc/h\n"%(sm)+\
				'%12s\t%12s'%('k','amplitude ratio average'),delimiter='\t')

#			cross1=cross_angle.cross_cor_angle(cor_avg1,delta_k2_avg_ini,delta_k2_avg)
#			cross2=cross_angle.cross_cor_angle(cor_avg2,delta_k2_avg_ini,delta_k2_avg)
#			cross3=cross_angle.cross_cor_angle(cor_avg3,delta_k2_avg_ini,delta_k2_avg)

#			file=open("cross1_rsd_%i_run%i_z%i_sm%i_%i.txt"%(n,run_n,redshift,smooth_scale,i),"w")
#			for k in cross1:
#				file.write(str(k)+'\n')
#			file.close()

#			file=open("cross2_rsd_%i_run%i_z%i_sm%i_%i.txt"%(n,run_n,redshift,smooth_scale,i),"w")
#			for k in cross2:
#				file.write(str(k)+'\n')
#			file.close()
	
#			file=open("cross3_rsd_%i_run%i_z%i_sm%i_%i.txt"%(n,run_n,redshift,smooth_scale,i),"w")
#			for k in cross3:
#				file.write(str(k)+'\n')
#			file.close()

	#	print('start plot_power_spectrum')
	#	ps.plot_power_spectrum(central_k,delta_k2_avg,i)
#			dis.displacement_correlation(S_s_x,S_s_y,S_s_z,S_tru_x,S_tru_y,S_tru_z,box_size,n,associated_kx,associated_ky,associated_kz,n3,run_n,redshift,smooth_scale,i)
			stop = timeit.default_timer()
			print("Iteration %d : Time : %f"%(i,stop-start))
#		np.save("displacement_est_x_rsd_run%i_z%i_sm%i.dat.npy"%(run_n,redshift,smooth_scale),S_s_x)
#		np.save("displacement_est_y_rsd_run%i_z%i_sm%i.dat.npy"%(run_n,redshift,smooth_scale),S_s_y)
#		np.save("displacement_est_z_rsd_run%i_z%i_sm%i.dat.npy"%(run_n,redshift,smooth_scale),S_s_z)
			delta_k_current=np.array(np.fft.rfftn(delta_L_current))*(L/N)**3. #  fftn for running next iter
#	filename='/home/xc298/palmer_scratch/scratch60/reconstruction/output/'\
#		+dirname+'/large'+'/delta_r_hada_finaliter_space%s_ani%.1f_run%i_sm%i_z%i.dat'%(space,Cani,run_n,smooth_scale,redshift) # for fiducial
	filename='/home/xc298/project/output/'\
		+dirname+'/large'+'/delta_r_hada_finaliter_space%s_ani%.1f_run%i_sm%i_z%i.dat'%(space,Cani,run_n,smooth_scale,redshift)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.save(filename,delta_L_current)
