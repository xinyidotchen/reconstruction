def standard(delta_k_bf_rec,N,L,redshift,delta_k_ini,mono_pk_ini,simulation_identifier,recon_params,space,mass_lim_low,mass_lim_up,Cani,pkxibins_spec):
	import numpy as np
	import standard_steps_general as steps
	import power_spectrum_mono_general as ps
	import power_spectrum_multi_general_pkxi as psm
	import gk_rk_general as cos
	import os
	import rk_avg_later_general as rk_avg_later
	import amplitude_recovery_general as amp
	import correlation_angle as cor_angle


	associated_kx=np.fft.fftfreq(N,d=1)*2*np.pi*N/L # kx vector. Python gives (0,1,...)/N. k=2*pi*q/L, so multiply it by 2*pi*N/L
	associated_ky=np.fft.fftfreq(N,d=1)*2*np.pi*N/L
	associated_kz=np.fft.fftfreq(N,d=1)[0:int(N/2+1)]*2*np.pi*N/L


	ksize = (len(associated_kx), len(associated_ky), len(associated_kz))
	        # Define k grid and associated k
	kxyz = np.meshgrid(associated_kx, associated_ky, associated_kz, indexing='ij')
	associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
#	associated_k[0,0,0] = 1.0
	associated_k[0,0,0]=1.0e-10
	
	smooth_scale=recon_params.smoothing
	Omega_m=recon_params.OmegaM(redshift)
	b=recon_params.bias
	f=recon_params.f

	mu=kxyz[2]/(associated_k*(1+1e-20))

	if space == 'real':
		delta_k_bf_rec=delta_k_bf_rec/b
	if space == 'z':
		delta_k_bf_rec=delta_k_bf_rec/(b+f*(mu**2))


#	print("delta_k_bf_rec x-axis mean=",np.mean(delta_k_bf_rec,axis=0))
#	print("delta_k_bf_rec y-axis mean=",np.mean(delta_k_bf_rec,axis=1))
#	print("delta_k_bf_rec z-axis mean=",np.mean(delta_k_bf_rec,axis=2))
#	print("delta_k_bf_rec x-axis std=",np.std(delta_k_bf_rec,axis=0))
#	print("delta_k_bf_rec y-axis std=",np.std(delta_k_bf_rec,axis=1))
#	print("delta_k_bf_rec z-axis std=",np.std(delta_k_bf_rec,axis=2))

	s_k_x,s_k_y,s_k_z=steps.zeldovich_dis(delta_k_bf_rec,kxyz,associated_k,redshift,smooth_scale,Omega_m,f,space,Cani)
#	s_r_x=np.fft.irfftn(s_k_x)/(L/N)**4.
#	s_r_y=np.fft.irfftn(s_k_y)/(L/N)**4.
#	s_r_z=np.fft.irfftn(s_k_z)/(L/N)**4.
#	print("s_r_x mean=",np.mean(s_r_x))
#	print("s_r_y mean=",np.mean(s_r_y))
#	print("s_r_z mean=",np.mean(s_r_z))
#
#	print("s_r_x std=",np.std(s_r_x))
#	print("s_r_y std=",np.std(s_r_y))
#	print("s_r_z std=",np.std(s_r_z))
#	
	print("starting displaced_original")

	category=simulation_identifier.category
	subcategory=simulation_identifier.subcategory
	run_n=simulation_identifier.sim_num
	dirname=simulation_identifier.get_standard_DirName(redshift,smooth_scale)
	header2_standard=simulation_identifier.get_standard_header(redshift,space,smooth_scale)
	if category=='Halos':
		dirname=simulation_identifier.get_standard_Halos_DirName(redshift,mass_lim_low,mass_lim_up,smooth_scale)
		header2_standard=simulation_identifier.get_standard_Halos_header(redshift,space,mass_lim_low,mass_lim_up,smooth_scale)

	if category=='Snapshots':	
		mass1=steps.displaced_original_Quijote_Snapshot(redshift,L,N,s_k_x,s_k_y,s_k_z,category,subcategory,run_n,space)
		print("ending displaced_original")
	
	if category=='Halos':
		mass1=steps.displaced_original_Quijote_Halos(redshift,L,N,s_k_x,s_k_y,s_k_z,category,subcategory,run_n,space,mass_lim_low,mass_lim_up)
		print("ending displaced_original")

	mass2=steps.shifted_density_field(L,N,s_k_x,s_k_y,s_k_z,f,space)


	delta_r_rec=steps.reconstructed_density_field(mass1,mass2,N)


#	filename='/home/xc298/palmer_scratch/scratch60/reconstruction/output/'\
#		+dirname+'/large'+'/delta_r_standard_space%s_ani%.1f.dat'%(space,Cani)
#	os.makedirs(os.path.dirname(filename), exist_ok=True)
#	np.save(filename,delta_r_rec)

#	delta_k_rec=np.fft.rfftn(delta_r_rec)*(L/N)**3
	delta_k_rec=np.array(np.fft.fftn(delta_r_rec))*(L/N)**3. # changed from rfftn to fftn to run analyses

	k_pk,r_cf=ps.power_spectrum(N,L,delta_k_rec,pkxibins_spec)
	mono_pk_rec=k_pk.T[2] # w/o window

	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/mono_pk_standard_space%s_ani%.1f.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, k_pk, fmt=['%.12f','%.12f','%.12f','%.12f'],\
		header="Monopole power spectrum "+header2_standard+'\n%12s\t%12s\t%12s\t%12s'%('k','P(k) w window','P(k) w/o window','k counts'),delimiter='\t')

	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/mono_cf_standard_space%s_ani%.1f.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, r_cf, fmt=['%.12f','%.12f','%.12f'],\
		header="Monopole correlation function "+header2_standard+'\n%12s\t%12s\t%12s'%('r','xi(r)','r counts'),delimiter='\t')

	

	k_pk0_pk2_numPl3,r_xi0_xi2=psm.power_spectrum(N,L,delta_r_rec,'standard',numPl=3,pkxibins_spec=pkxibins_spec)
	k_pk0_pk2_numPl5,r_xi0_xi2=psm.power_spectrum(N,L,delta_r_rec,'standard',numPl=5,pkxibins_spec=pkxibins_spec)

	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/multi_pk_standard_space%s_ani%.1f_numPl3.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, k_pk0_pk2_numPl3, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f','%.12f'],\
		header="Mono, quad power spectrum "+header2_standard+'\n%12s\t%12s\t%12s\t%12s\t%12s\t%12s'%('k','P0(k) w/ window','P2(k) w/ window','P0(k) w/o window','P2(k) w/o window','k counts'),delimiter='\t')

	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/multi_pk_standard_space%s_ani%.1f_numPl5.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, k_pk0_pk2_numPl5, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f','%.12f'],\
		header="Mono, quad power spectrum "+header2_standard+'\n%12s\t%12s\t%12s\t%12s\t%12s\t%12s'%('k','P0(k) w/ window','P2(k) w/ window','P0(k) w/o window','P2(k) w/o window','k counts'),delimiter='\t')


	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/multi_cf_standard_space%s_ani%.1f.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, r_xi0_xi2, fmt=['%.12f','%.12f','%.12f','%.12f'],\
		header="Mono, quad correlation function "+header2_standard+'\n%12s\t%12s\t%12s\t%12s'%('r','xi0(r)','xi2(r)','r counts'),delimiter='\t')
	


	k_avg,gk_rec,rk_rec,numer=cos.cross_cor_coeff(delta_k_ini,mono_pk_ini,delta_k_rec,mono_pk_rec,L,N,pkxibins_spec)
	k_avg,rk_avg_later_rec=rk_avg_later.rk_avg_later_cross_cor_coeff(delta_k_ini,delta_k_rec,L,N,pkxibins_spec)
	k_gk_rk=np.array([k_avg,gk_rec,rk_rec,rk_avg_later_rec,numer]).T

	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/gk_rk_standard_space%s_ani%.1f.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, k_gk_rk, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f'],\
		header="G(k) and r(k) "+header2_standard+'\n%12s\t%12s\t%12s\t%12s\t%12s'%('k','g(k)','r(k)','r(k) avg later','numer'),delimiter='\t')


	k_avg,cor_avg1,cor_avg2,cor_avg3,pk_avg1,pk_avg2,pk_avg3=cor_angle.correlation_angle(delta_k_rec,delta_k_ini,L,N,pkxibins_spec)
	cor_w_angle=np.array([k_avg,cor_avg1,cor_avg2,cor_avg3,pk_avg1,pk_avg2,pk_avg3]).T
	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/cor_angle_numer_standard_space%s_ani%.1f.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, cor_w_angle, fmt=['%.12f','%.12f','%.12f','%.12f','%.12f','%.12f','%.12f'],\
		header="cross correlation, angle (numerator) "+header2_standard+\
		'%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s'%('k','mu<1/3 numer','1/3<mu<2/3 numer','mu>2/3 numer','mu<1/3 pk','1/2<mu<2/3 pk','mu>2/3 pk'),delimiter='\t')


	k_avg,amplitude_ratio=amp.amplitude_recovery(delta_k_ini,delta_k_rec,L,N,pkxibins_spec)
	k_amp=np.array([k_avg,amplitude_ratio]).T        

	filename='/home/xc298/project/reconstruction_project/output/'\
		+dirname+'/amplitude_ratio_space%s_ani%.1f.txt'%(space,Cani)
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	np.savetxt(filename, k_amp, fmt=['%.12f','%.12f'],\
		header="amplitude ratio "+header2_standard+'\n%12s\t%12s'%('k','amplitude ratio average'),delimiter='\t')
