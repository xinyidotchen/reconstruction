def read_in_files_Quijote(box_size,grid_size,category,subcategory,run_n,redshift,space,isnap,mass_lim_low,mass_lim_up):
    import readgadget
    import h5py
    import hdf5plugin
    import readfof
    import numpy as np
    import os
    import random
    
    if category=='Snapshots':

        z_dict_inv={0.0:4,0.5:3,1.0:2,2.0:1,3.0:0}
        snapnum=z_dict_inv[redshift]
#        snapshot=["/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i/snapdir_00%i/snap_00%i.%i"%(category,subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)]

#        snapshot=["/home/xc298/palmer_scratch/simulations/quijote/%s/%s/%i/snapdir_00%i/snap_00%i.%i"%(category,subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)]
        snapshot=["/home/xc298/project/simulations/quijote/%s/%s/%i/snapdir_00%i/snap_00%i.%i"%(category,subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)] # for LC_p, LC_m

#        snapshot=["/home/xc298/palmer_scratch/simulations/quijote/%s/%s/%i/snapdir_00%i/snap_00%i.%i"%(category,subcategory,run_n,snapnum,snapnum,ii) for ii in range(0,8)] # for fiducial
        ptype = 1 # CDM

    	# read position
#        random.seed(isnap)
        pos = readgadget.read_field(snapshot[isnap],"POS ", ptype)/1e3

        ##### for test of lower number density
#        tot_num_particles=np.shape(pos)[0]
#        frac=1000**2/1024**3/10 # 1e-4 Mpc-3
#        dense=1e-4
#        num_par=int(tot_num_particles*frac)
#        random_mask=random.sample(range(tot_num_particles),num_par)

#        xlist=pos[random_mask,0]
#        ylist=pos[random_mask,1]
#        zlist=pos[random_mask,2]




#        vol=box_size**3
	
#        num_par_vol=np.column_stack((tot_num_particles,num_par,vol))
#        with open('/home/xc298/project/reconstruction_project/output/Quijote/'\
#                +'%s/%s/%i/grid%i/z%.1f/bf_rec/'%(category,subcategory,run_n,grid_size,redshift)+'num_particles_%s_%1.0e.txt'%(space,dense),"a") as a_file:
#            np.savetxt(a_file,num_par_vol)

        #####

        xlist=pos[:,0]
        ylist=pos[:,1]
        zlist=pos[:,2]

        x=xlist/box_size*grid_size
        y=ylist/box_size*grid_size
        z=zlist/box_size*grid_size

        if space=='z':
            vel =readgadget.read_field(snapshot[isnap],"VEL ",ptype)
#            z_v= vel[random_mask,2]  # for test of lower number density
            z_v= vel[:,2]

            def E(z):
                Omega_m=0.3175 ### Quijote
                Omega_L=1.-Omega_m
                return np.sqrt(Omega_m*(1+redshift)**3+Omega_L)

            a=1./(1.+redshift)
            z=(zlist+z_v/a/100./E(z))/box_size*grid_size
        



    if category=='Halos':
    	
#        snapdir="/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i"%(category,subcategory,run_n)
        snapdir="/home/xc298/palmer_scratch/simulations/quijote/%s/%s/%i"%(category,subcategory,run_n)
        z_dict_inv={0.0:4,0.5:3,1.0:2,2.0:1,3.0:0}
        snapnum=z_dict_inv[redshift]
        FoF=readfof.FoF_catalog(snapdir, snapnum, long_ids=False, swap=False, SFR=False, read_IDs=False)

        print("length before cut=", len(FoF.GroupMass))

        cut_mask=(np.log10(FoF.GroupMass*1e10) < mass_lim_up) & (np.log10(FoF.GroupMass*1e10) >= mass_lim_low)

        xlist=FoF.GroupPos[cut_mask,0]/1e3
        ylist=FoF.GroupPos[cut_mask,1]/1e3
        zlist=FoF.GroupPos[cut_mask,2]/1e3

        num_halos=np.float(len(xlist))
        vol=box_size**3

        num_halos_vol=np.array([num_halos,vol])
        filename='/home/xc298/project/reconstruction_project/output/Quijote/'\
                +'%s/%s/%i/grid%i/z%.1f/bf_rec/%.2f-%.2f/'%(category,subcategory,run_n,grid_size,redshift,mass_lim_low,mass_lim_up)+'num_halos_%s.txt'%(space)
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        np.savetxt(filename,num_halos_vol,header='# halos in mass bin %.2f-%.2f'%(mass_lim_low,mass_lim_up)+'\n%s\t%s'%('num_halos','volume[Mpc^3]'),delimiter='\t')

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

    return x,y,z
