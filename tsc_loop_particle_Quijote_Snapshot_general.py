def tsc_loop_particle(box_size,grid_size,run_n,redshift,subcategory,space):
	import numpy as np
	import timeit
	import pygadget as pyg
	import readgadget
	import random
	
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

	mass=np.zeros(ngrid**3)	    

	

	z_dict_inv={0.0:4,0.5:3,1.0:2,2.0:1,3.0:0}
	snapnum=z_dict_inv[redshift]
	snapshot=["/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%i/snapdir_00%i/snap_00%i.%i"%(subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)]
	ptype = 1 # CDM
	
	# read positions
	for isnap in range(8):
	    pos = readgadget.read_field(snapshot[isnap],"POS ", ptype)/1e3
	    print("pos shape=",np.shape(pos))

	    tot_num_particles=np.shape(pos)[0]
	    frac=1000**2/1024**3 # 1e-3 Mpc-3
	    cut_mask=random.sample(range(tot_num_particles),int(tot_num_particles*frac))

	    xlist=pos[random_mask,0]
	    ylist=pos[random_mask,1]
	    zlist=pos[random_mask,2]

#	    xlist=pos[:,0]
#	    ylist=pos[:,1]
#	    zlist=pos[:,2]

	    x=xlist/box_size*grid_size
	    y=ylist/box_size*grid_size
	    z=zlist/box_size*grid_size

	    if space=='z':
	        vel =readgadget.read_field(snapshot[isnap],"VEL ",ptype)
	        z_v= vel[random_mask,2]
#	        z_v= vel[:,2]
	        a=1./(1.+redshift)
	        z=(zlist+z_v/a/100.)/box_size*grid_size
	    print("z mean=",np.mean(z))
	    print("z median=",np.median(z))
	    print("x mean=",np.mean(x))
	    print("x median=",np.median(x))

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
	    
	    
	    
	#     print(len(mass))
	    for dix in range(-1,2): # apply weights
	        for diy in range(-1,2):
	            for diz in range(-1,2):
	                ws=weights(dix,diy,diz)
	#                 print((ws))
	                ibin = bins(dix,diy,diz)
	#                 print((ibin))
	                mass+=np.bincount(ibin,weights=ws,minlength=ngrid**3)
	    
	
	return np.reshape(mass,(ngrid,ngrid,ngrid))







def tsc_loop_particle_ini(box_size,grid_size,run_n,subcategory):
	import numpy as np
	import timeit
	import pygadget as pyg

	
	ngrid=grid_size


	def periodic(ng, ix):
	    ww = np.nonzero(ix>=ng)
	    ix[ww]-=ng
	    return ix
	    

	def periodic2(ng,ixm):
#	    ww=np.nonzero((ixm+ng)==(ng-1))
	    ww = np.nonzero(ixm <= (-1))
	    ixm[ww]+=(ng)
	    return ixm

	def weights(dix,diy,diz):
	    weights=weightx[dix]*weighty[diy]*weightz[diz]
	    return weights
	        
	def bins(dix,diy,diz):
		bins=binx[dix]*ngrid**2 + biny[diy]*ngrid + binz[diz]
		return bins	

	mass=np.zeros(ngrid**3)	    
	# Set the filenames
#	fns = ["/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/Snapshots/%s/%i/ICs/ics.%i"%(subcategory,run_n,ii) for ii in range(0,512)] # 128 for Mnu_pp.
	fns = ["/home/xc298/palmer_scratch/simulations/quijote/Snapshots/%s/%i/ICs/ics.%i"%(subcategory,run_n,ii) for ii in range(0,512)]
	# Define lists
	
	for ifn in fns:
	    start=timeit.default_timer() 
	    xlist = []
	    ylist = []
	    zlist = []

	    sim = pyg.Simulation(ifn)
	    block = sim.read_block("pos","halo")

	    xlist.append(np.array(block.x))
	    ylist.append(np.array(block.y))
	    zlist.append(np.array(block.z))

	    x = np.concatenate(xlist)
	    y = np.concatenate(ylist)
	    z = np.concatenate(zlist)

	    x=x/1e3/box_size*grid_size
	    y=y/1e3/box_size*grid_size
	    z=z/1e3/box_size*grid_size

	# time it

	    #mass=np.zeros(ngrid**3) # mass grid
	    ix = np.int32(x+0.5)
	    iy = np.int32(y+0.5)  
	    iz = np.int32(z+0.5)

	    
	    ixp=np.int32(x+0.5)+1
	    ixm=np.int32(x+0.5)-1
	#     print(ixm)
	    ixp_weight=(0.5*(1.5-abs(x-(ix+1)))**2.) # weight for the second nearst grid points
	    ixm_weight=(0.5*(1.5-abs(x-(ix-1)))**2.)
	    ix_weight=(0.75-(abs(x-ix))**2.) # weight for the nearest grid point

	    iyp=np.int32(y+0.5)+1
	    iym=np.int32(y+0.5)-1

	    iyp_weight=(0.5*(1.5-abs(y-(iy+1)))**2.) # weight for the second nearst grid points
	    iym_weight=(0.5*(1.5-abs(y-(iy-1)))**2.)    
	    iy_weight=(0.75-(abs(y-iy))**2.) # weight for the nearest grid point

	    izp=np.int32(z+0.5)+1
	    izm=np.int32(z+0.5)-1

	    izp_weight=(0.5*(1.5-abs(z-(iz+1)))**2.) # weight for the second nearst grid points
	    izm_weight=(0.5*(1.5-abs(z-(iz-1)))**2.)
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
	    
	    weightx=[ixm_weight,ix_weight,ixp_weight]
	    weighty=[iym_weight,iy_weight,iyp_weight]
	    weightz=[izm_weight,iz_weight,izp_weight]
	    
	    binx=[ixm,ix,ixp]
	    biny=[iym,iy,iyp]
	    binz=[izm,iz,izp]
	    
	    
	    
	#     print(len(mass))
	    for dix in range(-1,2): # apply weights
	        for diy in range(-1,2):
	            for diz in range(-1,2):
	                ws=weights(dix,diy,diz)
	#                 print((ws))
	                ibin = bins(dix,diy,diz)
	#                 print((ibin))
	                mass+=np.bincount(ibin,weights=ws,minlength=ngrid**3)
	    
	    
	#     print(mass)
	#     print(len(mass))
	    
	    
	    stop = timeit.default_timer()

	    print(stop-start)
	
	return np.reshape(mass,(ngrid,ngrid,ngrid))






