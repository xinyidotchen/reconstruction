import numpy as np





def power_spectrum(grid_size,box_size,delta_k,pkxibins_spec):
    import timeit
    start=timeit.default_timer() # time it
    N=grid_size
    L=box_size

    kbins_min=pkxibins_spec.kbins_min
    kbins_max_excl=pkxibins_spec.kbins_max_excl
    kbins_width=pkxibins_spec.kbins_width
    rbins_min=pkxibins_spec.rbins_min
    rbins_max_excl=pkxibins_spec.rbins_max_excl
    rbins_width=pkxibins_spec.rbins_width

    associated_kx=np.fft.fftfreq(N,d=1)*2*np.pi*N/L # kx vector. Python gives (0,1,...)/N. k=2*pi*q/L, so multiply it by 2*pi*N/L
    associated_ky=np.fft.fftfreq(N,d=1)*2*np.pi*N/L
#    associated_kz=np.fft.fftfreq(N,d=1)[0:int(N/2+1)]*2*np.pi*N/L # kz vector, reduced to half of the number of elements of cooridinates because of complex number
    associated_kz=np.fft.fftfreq(N,d=1)*2*np.pi*N/L # to use fft instead of rfft

    ksize = (len(associated_kx), len(associated_ky), len(associated_kz))
    kxyz = np.meshgrid(associated_kx, associated_ky, associated_kz, indexing='ij')
    associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
    associated_k[0,0,0] = 1.0e-10

    delta_k2_element=np.real(delta_k*np.conj(delta_k))
    correlation_element=np.real(np.fft.ifftn(delta_k2_element/L**3.))/(L/N)**3. # changed from irfftn to ifftn # changed back # changed to ifftn + np.real

#    kbins=np.logspace(np.log10(0.01),np.log10(1),100) # log k bins for plot
#    kbins=np.linspace(0.01,1,100) # bin size 0.01 # Quijote bins
#    kbins=np.arange(0.,1.,0.005) # DESI bins
    kbins=np.arange(kbins_min,kbins_max_excl,kbins_width)
    
    #Wk2_TSC=np.zeros(ksize,dtype=np.complex128)
    k_N=np.pi/L*N # Nyquist frequency
#    Wk2_TSC=(np.sinc(np.pi/k_N*kxyz[0]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[1]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[2]/2./np.pi))**4. # CIC window function squared (power=2 for CIC)
    Wk2_TSC=(np.sinc(np.pi/k_N*kxyz[0]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[1]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[2]/2./np.pi))**6. # TSC window function squared (power=3 for TSC)
#    Wk2_TSC=np.ones(ksize,dtype=np.float)

    associated_k=associated_k.ravel()
    delta_k2_element=delta_k2_element.ravel()
    Wk2_TSC=Wk2_TSC.ravel()
    count_ind=np.digitize(associated_k,bins=kbins)
    ws = associated_k*0.0+1.0
    
    count=np.bincount(count_ind,weights=ws,minlength=len(kbins))
    print("count=",count)
    #average k
    k_counts=np.bincount(count_ind,weights=ws*associated_k,minlength=len(kbins))
    k_avg_full=k_counts/count
    
    delta_k2_element=delta_k2_element/Wk2_TSC
    delta_k2=np.bincount(count_ind,weights=ws*delta_k2_element,minlength=len(kbins))

    delta_k2_avg=delta_k2/count/L**3.
    
    central_k=[(kbins[x]+kbins[x+1])/2. for x in range(len(kbins)-1)]
    
    #associated_x=np.fft.fftfreq(delta_r.shape[0],d=1)*2*np.pi*N/L
    
############################################################################
########### a version of Pk w/o window correction ##########################

    Wk2_none=np.ones(ksize,dtype=np.float).ravel()
    delta_k2_element=np.real(delta_k*np.conj(delta_k)).ravel()
    delta_k2_element=delta_k2_element/Wk2_none
    delta_k2=np.bincount(count_ind,weights=ws*delta_k2_element,minlength=len(kbins))
    delta_k2_avg_Wk2_none=delta_k2/count/L**3.

############################################################################


    correlation_element=correlation_element.ravel()
    associated_x=np.fft.fftfreq(N,d=1)*L
    associated_y=np.fft.fftfreq(N,d=1)*L
    associated_z=np.fft.fftfreq(N,d=1)*L


    rsize = (len(associated_x), len(associated_y), len(associated_z))     
    rxyz = np.meshgrid(associated_x, associated_y, associated_z, indexing='ij')
    associated_r=np.sqrt(rxyz[0]**2 + rxyz[1]**2 + rxyz[2]**2)
    associated_r=associated_r.ravel()
    #rbins=np.logspace(np.log10(10),np.log10(300),25)
#    rbins=np.arange(0,250,5) # DESI bin width
#    rbins=np.linspace(10,250,150) # quijote bin width
    rbins=np.arange(rbins_min,rbins_max_excl,rbins_width)

    count_r_ind=np.digitize(associated_r,bins=rbins)
    ws_r = associated_r*0.0+1.0
    count_r=np.bincount(count_r_ind,weights=ws_r,minlength=len(rbins)+1)

    #avg r
    r_counts=np.bincount(count_r_ind,weights=ws_r*associated_r,minlength=len(rbins))
    r_avg=r_counts/count_r

    correlation=np.bincount(count_r_ind,weights=ws_r*correlation_element,minlength=len(rbins)+1)
    correlation_avg=correlation/count_r
    central_r=[(rbins[x]+rbins[x+1])/2. for x in range(len(rbins)-1)]   
 


    stop=timeit.default_timer() # time it
    print(stop-start)
#    return k_avg_full,r_avg,delta_k2_avg,correlation_avg,count,count_r
    savearr_pk=np.array([k_avg_full,delta_k2_avg,delta_k2_avg_Wk2_none,count]).T
    savearr_xi=np.array([r_avg,correlation_avg,count_r]).T
    return savearr_pk,savearr_xi



