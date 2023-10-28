def tsc_loop_particle(grid_size,mass,x,y,z):
	import numpy as np
	import timeit

	
	ngrid=grid_size

	#periodic boundry condition for ix,ixp>=grid size: grid size -> 0, grid size+1 -> 1
	def periodic(ng, ix):
	    ww = np.nonzero(ix>=ng)
	    ix[ww]-=ng
	    return ix
	    
	#periodic boundry condition for ixm=-1: -1 -> grid size-1
	def periodic2(ng,ixm):
#	    ww=np.nonzero((ixm+ng)==(ng-1))
	    ww=np.nonzero((ixm)<=(-1))
	    ixm[ww]+=(ng)
	    return ixm

	def weights(dix,diy,diz):
	    weights=weightx[dix]*weighty[diy]*weightz[diz]
	    return weights
	#bins in 1D        
	def bins(dix,diy,diz):
		bins=binx[dix]*ngrid**2 + biny[diy]*ngrid + binz[diz]
		return bins	


	

	ix = np.int32(x+0.5)
	iy = np.int32(y+0.5)  
	iz = np.int32(z+0.5)

	    
	ixp=ix+1
	ixm=ix-1

	ixp_weight=(0.5*(1.5-abs(x-ixp))**2.) # weight for the second nearst grid points
	ixm_weight=(0.5*(1.5-abs(x-ixm))**2.)
	ix_weight=(0.75-(abs(x-ix))**2.) # weight for the nearest grid point

	iyp=iy+1
	iym=iy-1

	iyp_weight=(0.5*(1.5-abs(y-iyp))**2.) # weight for the second nearst grid points
	iym_weight=(0.5*(1.5-abs(y-iym))**2.)    
	iy_weight=(0.75-(abs(y-iy))**2.) # weight for the nearest grid point

	izp=iz+1
	izm=iz-1

	izp_weight=(0.5*(1.5-abs(z-izp))**2.) # weight for the second nearst grid points
	izm_weight=(0.5*(1.5-abs(z-izm))**2.)
	iz_weight=(0.75-(abs(z-iz))**2.) # weight for the nearest grid point

	ix = periodic(ngrid,ix)
	iy = periodic(ngrid,iy)
	iz = periodic(ngrid,iz)
	    
	ixm=periodic2(ngrid,ixm)
	ixp=periodic(ngrid,ixp)

	iym=periodic2(ngrid,iym)
	iyp=periodic(ngrid,iyp)
	    
	izm=periodic2(ngrid,izm)
	izp=periodic(ngrid,izp)

	ix=periodic2(ngrid,ix)
	iy=periodic2(ngrid,iy)
	iz=periodic2(ngrid,iz)

	ixm=periodic(ngrid,ixm)
	ixp=periodic2(ngrid,ixp)
	iym=periodic(ngrid,iym)
	iyp=periodic2(ngrid,iyp)
	izm=periodic(ngrid,izm)
	izp=periodic2(ngrid,izp)
	    

	weightx=[ixm_weight,ix_weight,ixp_weight]
	weighty=[iym_weight,iy_weight,iyp_weight]
	weightz=[izm_weight,iz_weight,izp_weight]
	    
	binx=[ixm,ix,ixp]
	biny=[iym,iy,iyp]
	binz=[izm,iz,izp]
	    
	    
	    
	for dix in range(-1,2): # apply weights
	    for diy in range(-1,2):
	        for diz in range(-1,2):
	            ws=weights(dix,diy,diz)
	#                 print((ws))
	            ibin = bins(dix,diy,diz)
	#                 print((ibin))
	            mass+=np.bincount(ibin,weights=ws,minlength=ngrid**3)
	    
	
	return np.reshape(mass,(ngrid,ngrid,ngrid))












