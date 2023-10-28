def do_tsc_Quijote(L,N,category,subcategory,run_n,redshift,space,mass_lim_low,mass_lim_up):
    import tsc_general as tsc
    import read_in_files_Quijote as read_in
    import numpy as np

    mass=np.zeros(N**3)
    if category == 'Snapshots':
        for isnap in range(8):
            x,y,z=read_in.read_in_files_Quijote(L,N,category,subcategory,run_n,redshift,space,isnap,mass_lim_low,mass_lim_up)
            mass=tsc.tsc_loop_particle(N,mass,x,y,z)
            mass=np.reshape(mass,(N**3))
#    mass=np.reshape(mass,(N,N,N))
    isnap=0
    if category == 'Halos':
        x,y,z=read_in.read_in_files_Quijote(L,N,category,subcategory,run_n,redshift,space,isnap,mass_lim_low,mass_lim_up)
        mass=tsc.tsc_loop_particle(N,mass,x,y,z)


    return np.reshape(mass,(N,N,N))
