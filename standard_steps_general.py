def zeldovich_dis(delta_k,kxyz,associated_k,redshift,smooth_scale,Omega_m,f,space,Cani):
	import smoothing_ani as sm
	import numpy as np
#	G_element=sm.G(associated_k,smooth_scale,Cani)
	G_element=sm.G(kxyz,smooth_scale,Cani)
	s_k_x=-1j*kxyz[0]*delta_k*G_element/associated_k**2
	s_k_y=-1j*kxyz[1]*delta_k*G_element/associated_k**2
	s_k_z=-1j*kxyz[2]*delta_k*G_element/associated_k**2

#	print("s_k_x mean=",np.mean(s_k_x))
#	print("s_k_y mean=",np.mean(s_k_y))
#	print("s_k_z mean=",np.mean(s_k_z))

#	print("s_k_x std=",np.std(s_k_x))
#	print("s_k_y std=",np.std(s_k_y))
#	print("s_k_z std=",np.std(s_k_z))
	if space=='z':
		s_k_z=s_k_z*(1+f)
	return s_k_x,s_k_y,s_k_z




def displaced_original_Quijote_Snapshot(redshift,box_size,grid_size,s_k_x,s_k_y,s_k_z,category,subcategory,run_n,space):
	from scipy.interpolate import RegularGridInterpolator
	from scipy.interpolate import interpn
	import numpy as np
	import pygadget as pyg
	import tsc_new as tsc
	import readfof
	import readgadget	
	import random
	
	
	N=grid_size
	L=box_size
	ngrid=grid_size
	
	s_r_x=np.fft.irfftn(s_k_x)/(L/N)**4.
	s_r_y=np.fft.irfftn(s_k_y)/(L/N)**4.
	s_r_z=np.fft.irfftn(s_k_z)/(L/N)**4.

	xgrid=np.linspace(0,511,512)	
	ygrid=np.linspace(0,511,512)
	zgrid=np.linspace(0,511,512)
	print(s_r_x.shape)
	inter_fun_x= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_x,bounds_error=False,fill_value=None)
	inter_fun_y= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_y,bounds_error=False,fill_value=None)
	inter_fun_z= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_z,bounds_error=False,fill_value=None)
	
	#print("inter function constructed")
	
	def periodic(ng, ix):
		ww = np.nonzero(ix>=ng)
		ix[ww]-=ng
		return ix

        #periodic boundry condition for ixm=-1: -1 -> grid size-1
	def periodic2(ng,ixm):
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
	def periodic3(ng,x0_chunck):
		ww1=np.nonzero(x0_chunck>ng/2.)
		x0_chunck[ww1]=x0_chunck[ww1]-ng
		ww2=np.nonzero(x0_chunck<-ng/2.)
		x0_chunck[ww2]=x0_chunck[ww2]+ng
		return x0_chunck

	mass1=np.zeros(N**3)

	
#	snapdir="/gpfs/loomis/scratch60/padmanabhan/xc298/reconstruction/code+input/large_box_copy10/large_box_copy9/Quijote_data/SD/10001"
#	snapnum=4
#	z_dict={4:0.0,3:0.5,2:1.0,0:3.0}
#	redshift=z_dict[snapnum]
#	FoF=readfof.FoF_catalog(snapdir, snapnum, long_ids=False,swap=False, SFR=False, read_IDs=False)
##	pos_h=FoF.GroupPos/1e3
	
##	xlist=pos_h[:,0]
##	ylist=pos_h[:,1]
##	zlist=pos_h[:,2]

#	cut_mask=(np.log10(FoF.GroupMass*1e10) < 13.75) & (np.log10(FoF.GroupMass*1e10) > 13.25)
#	xlist=FoF.GroupPos[cut_mask,0]/1e3	
#	ylist=FoF.GroupPos[cut_mask,1]/1e3
#	zlist=FoF.GroupPos[cut_mask,2]/1e3


##	vel_h=FoF.GroupVel*(1.0+redshift)

##	z_v=vel_h[:,2]
#	z_v=FoF.GroupVel[cut_mask,2]*(1.0+redshift)	






	z_dict_inv={0.0:4,0.5:3,1.0:2,2.0:1,3.0:0}
	snapnum=z_dict_inv[redshift]
#	snapshot=["/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i/snapdir_00%i/snap_00%i.%i"%(category,subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)]
	snapshot=["/home/xc298/palmer_scratch/simulations/quijote/%s/%s/%i/snapdir_00%i/snap_00%i.%i"%(category,subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)]
	ptype=1
	for isnap in range(8):
#		random.seed(isnap)
		pos  = readgadget.read_field(snapshot[isnap],"POS ", ptype)/1e3
		
		##### for test of lower number density
#		tot_num_particles=np.shape(pos)[0]
#		frac=1000**2/1024**3/10 # 1e-4 Mpc-3
#		num_par=int(tot_num_particles*frac)
#		random_mask=random.sample(range(tot_num_particles),num_par)

#		xlist=pos[random_mask,0]
#		ylist=pos[random_mask,1]
#		zlist=pos[random_mask,2]
		#####
		
		xlist= pos[:,0]
		ylist= pos[:,1]
		zlist= pos[:,2]
	
		x=xlist/box_size*grid_size
		y=ylist/box_size*grid_size
		z=zlist/box_size*grid_size

		if space == 'z':

			vel= readgadget.read_field(snapshot[isnap],"VEL ",ptype)
#			z_v= vel[random_mask,2]
			z_v=vel[:,2]
			a=1./(1.+redshift)
			def E(z):
				Omega_m=0.3175 ### Quijote
				Omega_L=1.-Omega_m
				return np.sqrt(Omega_m*(1+redshift)**3+Omega_L)
			z=(zlist+z_v/a/100./E(z))/box_size*grid_size	


		pts=np.column_stack((x,y,z))
		print("length of pts",len(pts))


		displaced_x=inter_fun_x(pts)
		displaced_y=inter_fun_y(pts)
		displaced_z=inter_fun_z(pts)
		
		print("displaced_x mean=", np.mean(displaced_x))
		print("displaced_y mean=", np.mean(displaced_y))
		print("displaced_z mean=", np.mean(displaced_z))

		print("displaced_x std=", np.std(displaced_x))
		print("displaced_y std=", np.std(displaced_y))
		print("displaced_z std=", np.std(displaced_z))
		
		

		x=x+displaced_x
		y=y+displaced_y
		z=z+displaced_z
	
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

					ibin = bins(dix,diy,diz)
					print("min ws=",np.min(ws))
					print("min ibin=",np.min(ibin))
					mass1+=np.bincount(ibin,weights=ws,minlength=ngrid**3)


	return np.reshape(mass1,(ngrid,ngrid,ngrid))






def displaced_original_Quijote_Halos(redshift,box_size,grid_size,s_k_x,s_k_y,s_k_z,category,subcategory,run_n,space,mass_lim_low,mass_lim_up):
        from scipy.interpolate import RegularGridInterpolator
        from scipy.interpolate import interpn
        import numpy as np
        import pygadget as pyg
        import tsc_new as tsc
        import readfof
        import readgadget
	
        N=grid_size
        L=box_size
        ngrid=grid_size

        s_r_x=np.fft.irfftn(s_k_x)/(L/N)**4.
        s_r_y=np.fft.irfftn(s_k_y)/(L/N)**4.
        s_r_z=np.fft.irfftn(s_k_z)/(L/N)**4.

        xgrid=np.linspace(0,511,512)
        ygrid=np.linspace(0,511,512)
        zgrid=np.linspace(0,511,512)
        print(s_r_x.shape)
        inter_fun_x= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_x,bounds_error=False,fill_value=None)
        inter_fun_y= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_y,bounds_error=False,fill_value=None)
        inter_fun_z= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_z,bounds_error=False,fill_value=None)

        #print("inter function constructed")

        def periodic(ng, ix):
                ww = np.nonzero(ix>=ng)
                ix[ww]-=ng
                return ix

        #periodic boundry condition for ixm=-1: -1 -> grid size-1
        def periodic2(ng,ixm):
                ww=np.nonzero((ixm)<=(-1))
                ixm[ww]+=(ng)
                return ixm

        def weights(dix,diy,diz):
                weights=weightx[dix]*weighty[diy]*weightz[diz]
                return weights
        def bins(dix,diy,diz):
                bins=binx[dix]*ngrid**2 + biny[diy]*ngrid + binz[diz]
                return bins
        def periodic3(ng,x0_chunck):
                ww1=np.nonzero(x0_chunck>ng/2.)
                x0_chunck[ww1]=x0_chunck[ww1]-ng
                ww2=np.nonzero(x0_chunck<-ng/2.)
                x0_chunck[ww2]=x0_chunck[ww2]+ng
                return x0_chunck

        mass1=np.zeros(N**3)

        snapdir="/home/xc298/palmer_scratch/simulations/quijote/%s/%s/%i"%(category,subcategory,run_n)
#        snapdir="/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i"%(category,subcategory,run_n)
        z_dict_inv={0.0:4,0.5:3,1.0:2,2.0:1,3.0:0}
        snapnum=z_dict_inv[redshift]
        FoF=readfof.FoF_catalog(snapdir, snapnum, long_ids=False,swap=False, SFR=False, read_IDs=False)


        cut_mask=(np.log10(FoF.GroupMass*1e10) < mass_lim_up) & (np.log10(FoF.GroupMass*1e10) >= mass_lim_low)
        xlist=FoF.GroupPos[cut_mask,0]/1e3      
        ylist=FoF.GroupPos[cut_mask,1]/1e3
        zlist=FoF.GroupPos[cut_mask,2]/1e3




        x=xlist/box_size*grid_size
        y=ylist/box_size*grid_size
        z=zlist/box_size*grid_size


        if space=='z':
            z_v=FoF.GroupVel[cut_mask,2]*(1.0+redshift)
            a=1./(1.+redshift)
            def E(z):
                Omega_m=0.3175 ### Quijote
                Omega_L=1.-Omega_m
                return np.sqrt(Omega_m*(1+redshift)**3+Omega_L)

            z=(zlist+z_v/a/100./E(z))/box_size*grid_size

        pts=np.column_stack((x,y,z))
        print("length of pts",len(pts))


        displaced_x=inter_fun_x(pts)
        displaced_y=inter_fun_y(pts)
        displaced_z=inter_fun_z(pts)

        print("displaced_x mean=", np.mean(displaced_x))
        print("displaced_y mean=", np.mean(displaced_y))
        print("displaced_z mean=", np.mean(displaced_z))

        print("displaced_x std=", np.std(displaced_x))
        print("displaced_y std=", np.std(displaced_y))
        print("displaced_z std=", np.std(displaced_z))



        x=x+displaced_x
        y=y+displaced_y
        z=z+displaced_z

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

                                ibin = bins(dix,diy,diz)
                                print("min ws=",np.min(ws))
                                print("min ibin=",np.min(ibin))
                                mass1+=np.bincount(ibin,weights=ws,minlength=ngrid**3)
        return np.reshape(mass1,(ngrid,ngrid,ngrid))



def shifted_density_field(box_size,grid_size,s_k_x,s_k_y,s_k_z,f,space):
	from scipy.interpolate import RegularGridInterpolator
	import numpy as np
	import itertools
	import numpy.random as npr
	# generate particles
	npr.seed(0)	

	L=box_size
	N=grid_size
	ngrid=grid_size	

	if space=='z':
		s_k_z=s_k_z/(1+f)

	s_r_x=np.fft.irfftn(s_k_x)/(L/N)**4.
	s_r_y=np.fft.irfftn(s_k_y)/(L/N)**4.
	s_r_z=np.fft.irfftn(s_k_z)/(L/N)**4.
	
	xgrid=np.linspace(0,511,512)	
	ygrid=np.linspace(0,511,512)
	zgrid=np.linspace(0,511,512)

	inter_fun_x= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_x,bounds_error=False,fill_value=None)
	inter_fun_y= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_y,bounds_error=False,fill_value=None)
	inter_fun_z= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_z,bounds_error=False,fill_value=None)

	def periodic(ng, ix):
		ww = np.nonzero(ix>=ng)
		ix[ww]-=ng
		return ix

        #periodic boundry condition for ixm=-1: -1 -> grid size-1
	def periodic2(ng,ixm):
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
	def periodic3(ng,x0_chunck):
		ww1=np.nonzero(x0_chunck>ng/2.)
		x0_chunck[ww1]=x0_chunck[ww1]-ng
		ww2=np.nonzero(x0_chunck<-ng/2.)
		x0_chunck[ww2]=x0_chunck[ww2]+ng
		return x0_chunck

	mass2=np.zeros(N**3)

	for i in range(8): # for halos, continue 512^3 randoms.
		x=npr.uniform(0,N,N**3)
		y=npr.uniform(0,N,N**3)
		z=npr.uniform(0,N,N**3)


#shift by s
	#x0_chunck=[x0[i:i + 134217728] for i in range(0, len(dis_x), 134217728)]
	#y0_chunck=[y0[i:i + 134217728] for i in range(0, len(dis_y), 134217728)]
	#z0_chunck=[z0[i:i + 134217728] for i in range(0, len(dis_z), 134217728)]



		pts=np.column_stack((x,y,z))

		#displaced_x=periodic3(N,inter_fun_x(pts))
		#displaced_y=periodic3(N,inter_fun_y(pts))
		#displaced_z=periodic3(N,inter_fun_z(pts))
		displaced_x=inter_fun_x(pts)
		displaced_y=inter_fun_y(pts)
		displaced_z=inter_fun_z(pts)

		print("displaced_x mean=", np.mean(displaced_x))
		print("displaced_y mean=", np.mean(displaced_y))
		print("displaced_z mean=", np.mean(displaced_z))

		print("displaced_x std=", np.std(displaced_x))
		print("displaced_y std=", np.std(displaced_y))
		print("displaced_z std=", np.std(displaced_z))
		x=x+displaced_x
		y=y+displaced_y
		z=z+displaced_z

		#shifted_x=x+inter_fun_x(pts)
		#shifted_y=y+inter_fun_y(pts)
		#shifted_z=z+inter_fun_z(pts)
		
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

					ibin = bins(dix,diy,diz)

					mass2+=np.bincount(ibin,weights=ws,minlength=ngrid**3)
		#return shifted_x,shifted_y,shifted_z
		return np.reshape(mass2,(ngrid,ngrid,ngrid))

#mass1=tsc.tsc_new(ngrid,displaced_x,displaced_y,displaced_z)
#mass2=tsc.tsc_new(ngrid,shifted_x,shifted_y,shifted_z)

def reconstructed_density_field(mass1,mass2,grid_size):
	import numpy as np
	grid_vol=grid_size**3 # grid volume
	avg_rho1=np.sum(mass1)/grid_vol # average density rho
	delta_d=(mass1/1.-avg_rho1)/avg_rho1 # delta_r=(rho-average rho)/average rho


	avg_rho2=np.sum(mass2)/grid_vol
	delta_s=(mass2/1.-avg_rho2)/avg_rho2

	delta_r=delta_d-delta_s

	return delta_r




def displaced_original_desi(redshift,box_size,grid_size,s_k_x,s_k_y,s_k_z,stage,SIM,tracer,catalog,space):
        from scipy.interpolate import RegularGridInterpolator
        from scipy.interpolate import interpn
        import numpy as np
        import pygadget as pyg
        import tsc_new as tsc
        import readfof
        import readgadget
	
        N=grid_size
        L=box_size
        ngrid=grid_size

        s_r_x=np.fft.irfftn(s_k_x)/(L/N)**4.
        s_r_y=np.fft.irfftn(s_k_y)/(L/N)**4.
        s_r_z=np.fft.irfftn(s_k_z)/(L/N)**4.

        xgrid=np.linspace(0,511,512)
        ygrid=np.linspace(0,511,512)
        zgrid=np.linspace(0,511,512)
        print(s_r_x.shape)
        inter_fun_x= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_x,bounds_error=False,fill_value=None)
        inter_fun_y= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_y,bounds_error=False,fill_value=None)
        inter_fun_z= RegularGridInterpolator((xgrid, ygrid, zgrid), s_r_z,bounds_error=False,fill_value=None)

        #print("inter function constructed")

        def periodic(ng, ix):
                ww = np.nonzero(ix>=ng)
                ix[ww]-=ng
                return ix

        #periodic boundry condition for ixm=-1: -1 -> grid size-1
        def periodic2(ng,ixm):
                ww=np.nonzero((ixm)<=(-1))
                ixm[ww]+=(ng)
                return ixm

        def weights(dix,diy,diz):
                weights=weightx[dix]*weighty[diy]*weightz[diz]
                return weights
        def bins(dix,diy,diz):
                bins=binx[dix]*ngrid**2 + biny[diy]*ngrid + binz[diz]
                return bins
        def periodic3(ng,x0_chunck):
                ww1=np.nonzero(x0_chunck>ng/2.)
                x0_chunck[ww1]=x0_chunck[ww1]-ng
                ww2=np.nonzero(x0_chunck<-ng/2.)
                x0_chunck[ww2]=x0_chunck[ww2]+ng
                return x0_chunck

        mass1=np.zeros(N**3)

######### Temporary, for MCI ##########
        xlist = np.genfromtxt("/home/xc298/scratch60/reconstruction/code+input/folder_for_large_box_copy10_large_box_copy9/DESI_data/cosmosim_UNIT-BAO-RSD-challenge_1Gpc_ELG-single-tracer/UNIT_DESI_Shadab_HOD_snap97_ELG_v1_4col.txt",usecols=0,unpack='True')
        ylist = np.genfromtxt("/home/xc298/scratch60/reconstruction/code+input/folder_for_large_box_copy10_large_box_copy9/DESI_data/cosmosim_UNIT-BAO-RSD-challenge_1Gpc_ELG-single-tracer/UNIT_DESI_Shadab_HOD_snap97_ELG_v1_4col.txt",usecols=1,unpack='True')
        zlist = np.genfromtxt("/home/xc298/scratch60/reconstruction/code+input/folder_for_large_box_copy10_large_box_copy9/DESI_data/cosmosim_UNIT-BAO-RSD-challenge_1Gpc_ELG-single-tracer/UNIT_DESI_Shadab_HOD_snap97_ELG_v1_4col.txt",usecols=2,unpack='True')
#######################################


#       xlist=np.genfromtxt("/home/xc298/scratch60/reconstruction/DESI/data/%s/%s/%s/%s"%(stage,SIM,tracer,catalog),usecols=0,unpack='True')
#       ylist=np.genfromtxt("/home/xc298/scratch60/reconstruction/DESI/data/%s/%s/%s/%s"%(stage,SIM,tracer,catalog),usecols=1,unpack='True')
#       zlist=np.genfromtxt("/home/xc298/scratch60/reconstruction/DESI/data/%s/%s/%s/%s"%(stage,SIM,tracer,catalog),usecols=2,unpack='True')





        x=xlist/box_size*grid_size
        y=ylist/box_size*grid_size
        z=zlist/box_size*grid_size


        if space=='z':

######## Temporary, for MCI ############
                zlist = np.genfromtxt("/home/xc298/scratch60/reconstruction/code+input/folder_for_large_box_copy10_large_box_copy9/DESI_data/cosmosim_UNIT-BAO-RSD-challenge_1Gpc_ELG-single-tracer/UNIT_DESI_Shadab_HOD_snap97_ELG_v1_4col.txt",usecols=3,unpack='True')
########################################

                z=zlist/box_size*grid_size




        pts=np.column_stack((x,y,z))
        print("length of pts",len(pts))


        displaced_x=inter_fun_x(pts)
        displaced_y=inter_fun_y(pts)
        displaced_z=inter_fun_z(pts)

        print("displaced_x mean=", np.mean(displaced_x))
        print("displaced_y mean=", np.mean(displaced_y))
        print("displaced_z mean=", np.mean(displaced_z))

        print("displaced_x std=", np.std(displaced_x))
        print("displaced_y std=", np.std(displaced_y))
        print("displaced_z std=", np.std(displaced_z))



        x=x+displaced_x
        y=y+displaced_y
        z=z+displaced_z

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

                                ibin = bins(dix,diy,diz)
                                print("min ws=",np.min(ws))
                                print("min ibin=",np.min(ibin))
                                mass1+=np.bincount(ibin,weights=ws,minlength=ngrid**3)
        return np.reshape(mass1,(ngrid,ngrid,ngrid))







