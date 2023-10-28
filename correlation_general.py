def correlation(delta_k_2,delta_k_1,L,N,pkxibins_spec):
	import numpy as np

	kbins_min=pkxibins_spec.kbins_min
	kbins_max_excl=pkxibins_spec.kbins_max_excl
	kbins_width=pkxibins_spec.kbins_width


	correlation_delta_k_2_delta_k_1= np.real(np.conj(delta_k_2)*delta_k_1)

	associated_kx=np.fft.fftfreq(N,d=1)*2*np.pi*N/L # kx vector. Python gives (0,1,...)/N. k=2*pi*q/L, so multiply it by 2*pi*N/L
	associated_ky=np.fft.fftfreq(N,d=1)*2*np.pi*N/L
#	associated_kz=np.fft.fftfreq(N,d=1)[0:int(N/2+1)]*2*np.pi*N/L # kz vector, reduced to half of the number of elements of cooridinates because of complex number
	associated_kz=np.fft.fftfreq(N,d=1)*2*np.pi*N/L # to incorporate using rfftn now

	ksize = (len(associated_kx), len(associated_ky), len(associated_kz))
	kxyz = np.meshgrid(associated_kx, associated_ky, associated_kz, indexing='ij')
	associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
#	kbins=np.logspace(np.log10(0.01),np.log10(1),100) # log k bins for plot
#	kbins=np.linspace(0.01,1,100) # quijote k bins
#	kbins=np.arange(0.,1.,0.005) # DESI k bins
	kbins=np.arange(kbins_min,kbins_max_excl,kbins_width)

	associated_k=associated_k.ravel()
	correlation_delta_k_2_delta_k_1=correlation_delta_k_2_delta_k_1.ravel()
	count_ind=np.digitize(associated_k,bins=kbins)
	ws = associated_k*0.0+1.0
	count=np.bincount(count_ind,weights=ws,minlength=len(kbins))

	k_counts=np.bincount(count_ind,weights=ws*associated_k,minlength=len(kbins))
	k_avg_full=k_counts/count

	correlation_sum=np.bincount(count_ind,weights=ws*correlation_delta_k_2_delta_k_1,minlength=len(kbins))
	correlation_avg=correlation_sum/count/L**3.
	return k_avg_full,correlation_avg
