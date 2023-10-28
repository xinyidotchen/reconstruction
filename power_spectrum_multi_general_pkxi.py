import numpy as np
from npmisc import *

def power_spectrum(grid_size,box_size,delta_r,sim_identifier,numPl,pkxibins_spec):
    import timeit
    import scipy.special

    start=timeit.default_timer() # time it
    N=grid_size
    L=box_size

    kbins_min=pkxibins_spec.kbins_min
    kbins_max_excl=pkxibins_spec.kbins_max_excl
    kbins_width=pkxibins_spec.kbins_width
    rbins_min=pkxibins_spec.rbins_min
    rbins_max_excl=pkxibins_spec.rbins_max_excl
    rbins_width=pkxibins_spec.rbins_width


    ksize, kxyz, associated_k = make_k_grid(N, L)
    mu=kxyz[2]/associated_k

    delta_k=np.array(np.fft.fftn(delta_r))*(L/N)**3.  # use fftn instead of rfftn so that I don't need to adjust the average later.

    print("min mu=",np.min(mu))
    print("max mu=",np.max(mu))
    delta_k2_element=np.real(delta_k*np.conj(delta_k)).flatten()

#    kbins=np.arange(0.,1.,0.005)
    kbins=np.arange(kbins_min,kbins_max_excl,kbins_width)
    mubins=np.linspace(-1.,1.,200) # number of mu bins does not matter. Increased to 200 from 48 here so that central_mu more likely to represent the mean of mu in each bin.
    dmu = mu[1]-mu[0]
    mubins[0] -= 1.0e-10 # to have symmetric mu after the next step.
    mubins[-1] += 1.0e-10 # to include mu=1 in the bin with now right edge 1+1e-10, so that the last bin is not 1<=bin<inf which only has one element, 1.

    count_ind_k=np.digitize(associated_k,bins=kbins) # output bin index 0 has no element. -1 to shift index to start from 0. # using k bin starting at 0.01, then index 0 has elements, so removing "-1"
    count_ind_mu=np.digitize(mu,bins=mubins)-1 # -1 is to shift index to start from 0. Can do that because 0th bin is empty. By default, digitize add one bin before the first bin with index 0, whether or not the first bin covers the minimum value. If there is an element smaller than the left edge of the first bin, the 0th bin will have a non-zero count. If there is no such element, the 0th bin has zero count.

    count_ind=((count_ind_k)*(len(mubins)-1)+count_ind_mu).flatten() # len(mubins)-1 is to subtract the last bin, 1+1e-10<=mu bin<inf, which has no elements.
    max1dbins = (len(kbins)+1)*(len(mubins)-1) # len(mubins)-1 same reason as above. # len(kbins)+1 because now I use kbins (dk=0.01) starting from 0.01 instead of 0, so the 0th bin (k<0.01)has non-zero elements.
    kmushape = (len(kbins)+1, len(mubins)-1) # len(mubins)-1 same reason as above. # len(kbins)+1 same reason as above. To fix error msg from fiducial_HR_snap_main1_z=0_zspace_multi_kbin_0.01_15468293_0.out
    print(len(count_ind))
    ws=count_ind*0.0+1.0

    k_N=np.pi/L*N # Nyquist frequency
#    Wk2_TSC=((np.sinc(np.pi/k_N*kxyz[0]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[1]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[2]/2./np.pi))**4.).ravel() # CIC window function for CIC squared (power=2 for CIC)
    Wk2_TSC=((np.sinc(np.pi/k_N*kxyz[0]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[1]/2./np.pi)*np.sinc(np.pi/k_N*kxyz[2]/2./np.pi))**6.).ravel()  # window function for TSC 3D squared (power=3 for TSC)

    Wk2_none=np.ones(ksize,dtype=np.float).ravel()
    P0k_ww=[]
    P2k_ww=[]
    for Wk2 in [Wk2_TSC,Wk2_none]:
        delta_k2_element=np.real(delta_k*np.conj(delta_k)).flatten()
        delta_k2_element=delta_k2_element/Wk2 # apply window function.
    
        count=np.bincount(count_ind,weights=ws,minlength=max1dbins) # number of grid points in each k-mu bin
        delta_k2=np.bincount(count_ind,weights=delta_k2_element,minlength=max1dbins)
        mu_avg = np.bincount(count_ind,weights=mu.flatten(),minlength=max1dbins)
        k_avg = np.bincount(count_ind,weights=associated_k.flatten(),minlength=max1dbins)
    
        mask = (count < 0.001)
        count[mask] = 1
        mask2d = mask.reshape(kmushape)
        count2d = count.reshape(kmushape)
    
        delta_k2_avg_full=(delta_k2/count/L**3).reshape(kmushape)
        mu_avg = (mu_avg/count).reshape(kmushape)
        k_avg = k_avg.reshape(kmushape)
    
        P_k_mat = np.zeros((len(kbins),numPl),dtype='f8')
        k_mean = np.zeros(len(kbins), dtype='f8')
        k_counts = np.zeros(len(kbins), dtype='f8')
        nmuarr = np.zeros(len(kbins), dtype='f8') 
        for i in range(len(kbins)):
            mygood = ~mask2d[i,:].flatten()
            nmu = np.sum(mygood) # to get number of non-zero mu's.
            nmuarr[i] = nmu
            if (nmu/2 < numPl): # /2 because negative mu gives the same Legendre polynomial value as the positive mu. P(k,-mu)=P(k,+mu). So half the equations are redundent. There is always one warning for k=0? 
                print("WARNING : number of mu values < number of Legendre polynomials")
            central_mu = mu_avg[i,mygood].flatten()
            Legen_mat=np.array([scipy.special.legendre(i)(central_mu) for i in range(0,2*numPl,2)])
            # Orthogonalize wrt to the monopole, probably should also do this for the quadrupole and hexadecapole
            # for ll in range(1,numPl):
            #     Legen_mat[ll,:] -= np.mean(Legen_mat[ll,:])
            Legen_mat = Legen_mat.T
            myy =delta_k2_avg_full[i,mygood].flatten()
    
            mu_counts=count2d[i,mygood]
            weight =  np.array(np.sqrt(mu_counts)) # weigh the RHS by sqrt(n) where n is number of mu's in the k-mu bin
            weighted_y= myy*weight
            weighted_Legen_mat= Legen_mat*weight.reshape((len(weight),1))
    
            sol = np.linalg.lstsq(weighted_Legen_mat,weighted_y,rcond=None)[0]
            P_k_mat[i,:] = sol
            k_counts[i] = np.sum(count2d[i,mygood]) # sum counts over mu
            k_mean[i] = np.sum(k_avg[i,mygood])/k_counts[i] # average over mu
  
        P0k=np.array(P_k_mat).T[0]
        P2k=np.array(P_k_mat).T[1]
        P0k_ww.append(P0k)
        P2k_ww.append(P2k)

    P0k_ww=np.array(P0k_ww)
    P2k_ww=np.array(P2k_ww)

    savearr = np.array([k_mean, P0k_ww[0], P2k_ww[0],P0k_ww[1],P2k_ww[1], k_counts]).T

    stop=timeit.default_timer() # time it
    print(stop-start)



    ##### xil #####

    delta_k2_element=np.real(delta_k*np.conj(delta_k))
    correlation_element=np.real(np.fft.ifftn(delta_k2_element/L**3.))/(L/N)**3.

#    rbins=np.linspace(10,300,80)
#    rbins=np.arange(0,250,5) # DESI bin width
#    rbins=np.linspace(10,250,150) # quijote bin width
    rbins=np.arange(rbins_min,rbins_max_excl,rbins_width)

    mubins=np.linspace(-1.,1.,48) # number of mu bins does not matter. Increased to 200 from 48 here so that central_mu more likely to represent the mean of mu in each bin.
    mubins[0] -= 1.0e-10 # to have symmetric mu after the next step.
    mubins[-1] += 1.0e-10 # to include mu=1 in the bin with now right edge 1+1e-10, so that the last bin is not 1<=bin<inf which only has one element, 1.


    associated_r=associated_k.flatten()/(2*np.pi*N/L)*L # kx,ky,kz=2piN/L*fftfreq. x,y,z=L*fftfreq


    correlation_element=correlation_element.ravel()

    count_ind_r=np.digitize(associated_r,bins=rbins) # not -1 because rbins starts at 10 and there is r below 10, so there are element in the bin<10.
    count_ind_mu=np.digitize(mu.flatten(),bins=mubins)-1
    count_r_ind=(count_ind_r)*(len(mubins)-1)+(count_ind_mu)
    print(count_r_ind)

#       ws_r = associated_r*0.0+1.0
    ws_r=count_r_ind*0.0+1.0

    count_r_2d=np.bincount(count_r_ind,weights=ws_r,minlength=len(rbins)*(len(mubins)-1)) # 2d means r-mu grid. Here is 1d
    print("first 100 r bins along mu=", np.sum(count_r_2d.reshape(len(rbins)+1,(len(mubins)-1)),axis=1)[0:100]) # debug 200 mu bin
    print("count nonzero=",np.count_nonzero(count_r_2d))
    print("zero indices=",np.where(count_r_2d==0)[0:100])
    print(len(count_r_2d))
    print(np.max(count_r_2d))
    print(np.min(count_r_2d))
    print(np.mean(count_r_2d))
    del_rbins=int(np.ceil(np.max(np.where(count_r_2d==0))/(len(mubins)-1)))
    start_ind=int(del_rbins*(len(mubins)-1))
    count_r_2d=count_r_2d[start_ind:]

    print("count nonzero=",np.count_nonzero(count_r_2d))


    print("zero indices=",np.where(count_r_2d==0))
    print(len(count_r_2d))
    print(np.max(count_r_2d))
    print(np.min(count_r_2d))
    print(np.mean(count_r_2d))

    correlation=np.bincount(count_r_ind,weights=ws_r*correlation_element,minlength=len(rbins)*(len(mubins)-1))
    correlation=correlation[start_ind:]

    correlation_avg=correlation/count_r_2d
    correlation_avg=np.reshape(correlation_avg,(len(rbins)-del_rbins+1,(len(mubins)-1))) # +1 based on error message (slurm 6764886). why?

    r_avg=np.bincount(count_r_ind,weights=ws_r*associated_r,minlength=len(rbins)*(len(mubins)-1))
    r_avg=r_avg[start_ind:]

    r_avg=np.reshape(r_avg,(len(rbins)-del_rbins+1,(len(mubins)-1))) ## +1 based on error message (slurm 6764886). why? There are elements <10 and >250?    
    count_r_2d=np.reshape(count_r_2d,(len(rbins)-del_rbins+1,(len(mubins)-1))) ## +1 based on error message (slurm 6764886).
    r_mean=np.sum(r_avg,axis=1)/np.sum(count_r_2d,axis=1) # previously was taking average after average, which can be biased (if all points are in one r-mu bin, the average calculated that way woudl be smaller.)
    r_counts=np.sum(count_r_2d,axis=1)

    central_r=[(rbins[x]+rbins[x+1])/2. for x in range(len(rbins)-1)]

    central_mu=(mubins[1:]+mubins[0:-1])*0.5
    print("central_mu len",len(central_mu))
    print("mubins len",len(mubins))
    Legen2=(3.*np.array(central_mu)**2.-1.)/2.
    Legen0=np.ones(len(central_mu))
    d_mubin=mubins[1]-mubins[0]
    d_kbin=kbins[1]-kbins[0]


    xi2r=np.dot(correlation_avg, Legen2)*d_mubin/(np.dot(Legen2,Legen2)*d_mubin)
    xi0r=(1./(np.dot(Legen0,Legen0)*d_mubin))*np.dot(correlation_avg,Legen0)*d_mubin
    print("Legen orth coeff l=0:",np.dot(Legen0,Legen0)*d_mubin)
    print("Legen orth coeff l=2:",np.dot(Legen2,Legen2)*d_mubin)

    correlation_avg=np.reshape(correlation_avg,((len(mubins)-1)*((len(rbins)-del_rbins)+1))) ## +1 based on error message (slurm 6764886). why?
    print("r_mean len=",len(r_mean))
    print("xi0r len=",len(xi0r))
    print("xi2r len=",len(xi2r))
    print("r count len=",len(r_counts))
    savearr_xi = np.array([r_mean, xi0r,xi2r,r_counts]).T

    return savearr,savearr_xi
