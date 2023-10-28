from npmisc import *
def correlation_angle(delta_k_current,delta_k_ini,L,n,pkxibins_spec):
	import numpy as np


	kbins_min=pkxibins_spec.kbins_min
	kbins_max_excl=pkxibins_spec.kbins_max_excl
	kbins_width=pkxibins_spec.kbins_width
	
	cor_delta_k_nth_delta_k_ini= np.real(np.conj(delta_k_current)*delta_k_ini)



	        
#	ksize = (len(associated_kx), len(associated_ky), len(associated_kz))
#	kxyz = np.meshgrid(associated_kx, associated_ky, associated_kz, indexing='ij')
#	associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
#	associated_k[0,0,0]=1.0e-10
	        
#	mu=abs(kxyz[2]/associated_k)


	ksize, kxyz, associated_k = make_k_grid(n, L)
	mu=kxyz[2]/associated_k

	        
#	kbins=np.logspace(np.log10(0.01),np.log10(5),50) # log k bins for plot
#	kbins=np.linspace(0.01,1,100)        
#	kbins=np.arange(0.,1.,0.005)
	kbins=np.arange(kbins_min,kbins_max_excl,kbins_width)
	        
	associated_k=associated_k.ravel()
	mu=mu.ravel()
	cor_delta_k_nth_delta_k_ini=cor_delta_k_nth_delta_k_ini.ravel()

	        
	mask=~np.isnan(mu)
	mu=mu[mask]
	associated_k=associated_k[mask]
	cor_delta_k_nth_delta_k_ini=cor_delta_k_nth_delta_k_ini[mask]
	        
	mask1=np.abs(mu)<(1./3.)
	associated_k1=associated_k[mask1]
	cor_delta_k_nth_delta_k_ini1=cor_delta_k_nth_delta_k_ini[mask1]
	print("min mu=",np.min(mu))
	print("max mu=",np.max(mu))

	mask2=np.logical_and(np.abs(mu)>(1./3.),np.abs(mu)<(2./3))
	associated_k2=associated_k[mask2]
	cor_delta_k_nth_delta_k_ini2=cor_delta_k_nth_delta_k_ini[mask2]
	        
	mask3=np.abs(mu)>(2./3.)
	associated_k3=associated_k[mask3]
	cor_delta_k_nth_delta_k_ini3=cor_delta_k_nth_delta_k_ini[mask3]
	        
	        
	count_ind1=np.digitize(associated_k1,bins=kbins)
	ws1 = associated_k1*0.0+1.0
	count1=np.bincount(count_ind1,weights=ws1,minlength=len(kbins))
	cor_sum1=np.bincount(count_ind1,weights=ws1*cor_delta_k_nth_delta_k_ini1,minlength=len(kbins))
	cor_avg1=cor_sum1/count1/L**3.
	        
	        
	count_ind2=np.digitize(associated_k2,bins=kbins)
	ws2 = associated_k2*0.0+1.0
	count2=np.bincount(count_ind2,weights=ws2,minlength=len(kbins))
	cor_sum2=np.bincount(count_ind2,weights=ws2*cor_delta_k_nth_delta_k_ini2,minlength=len(kbins))
	cor_avg2=cor_sum2/count2/L**3.

	        
	count_ind3=np.digitize(associated_k3,bins=kbins)
	ws3 = associated_k3*0.0+1.0
	count3=np.bincount(count_ind3,weights=ws3,minlength=len(kbins))
	cor_sum3=np.bincount(count_ind3,weights=ws3*cor_delta_k_nth_delta_k_ini3,minlength=len(kbins))
	cor_avg3=cor_sum3/count3/L**3.

	count_ind=np.digitize(associated_k,bins=kbins)
	ws = associated_k*0.0+1.0
	count=np.bincount(count_ind,weights=ws,minlength=len(kbins))
	k_counts=np.bincount(count_ind,weights=ws*associated_k,minlength=len(kbins))
	k_avg=k_counts/count

	# file=open("cor_avg1_rsd_%i_run%i.txt"%(n,run_n),"w")
	# for i in cor_avg1:
	#         file.write(str(i)+'\n')
	# file.close()

	# file=open("cor_avg2_rsd_%i_run%i.txt"%(n,run_n),"w")
	# for i in cor_avg2:
	#         file.write(str(i)+'\n')
	# file.close()

	# file=open("cor_avg3_rsd_%i_run%i.txt"%(n,run_n),"w")
	# for i in cor_avg3:
	#         file.write(str(i)+'\n')
	# file.close()
	
	################## debug ########
#	print("count1=",count1)
#	print("count2=",count2)
#	print("count3=",count3)
#	print("tot count=",count)

	delta_k2_element=np.real(delta_k_current*np.conj(delta_k_current))
	delta_k2_element=delta_k2_element.ravel()
	delta_k2_element=delta_k2_element[mask]
	delta_k2_element1=delta_k2_element[mask1]
	delta_k2_element2=delta_k2_element[mask2]
	delta_k2_element3=delta_k2_element[mask3]

	pk_current_sum1=np.bincount(count_ind1,weights=ws1*delta_k2_element1,minlength=len(kbins))
	pk_avg1=pk_current_sum1/count1/L**3

	pk_current_sum2=np.bincount(count_ind2,weights=ws2*delta_k2_element2,minlength=len(kbins))
	pk_avg2=pk_current_sum2/count2/L**3

	pk_current_sum3=np.bincount(count_ind3,weights=ws3*delta_k2_element3,minlength=len(kbins))
	pk_avg3=pk_current_sum3/count3/L**3

	return k_avg,cor_avg1,cor_avg2,cor_avg3,pk_avg1,pk_avg2,pk_avg3
