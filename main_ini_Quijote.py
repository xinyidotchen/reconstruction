import numpy as np
import timeit
import tsc_loop_particle_Quijote_Snapshot_general as tsc
#import do_tsc_Quijote as tsc
import power_spectrum_mono_general as ps
import gk_rk_general as cos
import sys
from Quijote_names import Quijote_names
import os
from pkxibins_specifications import pkxibins_specifications

N=np.int(sys.argv[1]) # grid size
L=np.int(sys.argv[2]) # box size
category=np.str(sys.argv[3])
subcategory=np.str(sys.argv[4])
run_n=np.int(sys.argv[5]) # run number

kbins_min=np.float(sys.argv[6])
kbins_max_excl=np.float(sys.argv[7])
kbins_width=np.float(sys.argv[8])
rbins_min=np.float(sys.argv[9])
rbins_max_excl=np.float(sys.argv[10])
rbins_width=np.float(sys.argv[11])

n3=N**3

Quijote_identifier=Quijote_names(category,subcategory,run_n,N)

Quijote_ini_dirname=Quijote_identifier.get_ini_DirName()
sim_path=Quijote_identifier.getSimPath()


pkxibins_spec=pkxibins_specifications(kbins_min,kbins_max_excl,kbins_width,rbins_min,rbins_max_excl,rbins_width)




print("Starting ini TSC")
mass_ini=tsc.tsc_loop_particle_ini(L,N,run_n,subcategory)
print("Ending ini TSC")

grid_vol=n3 # grid volume
avg_rho=np.sum(mass_ini)/grid_vol # average density rho
delta_r_ini=(mass_ini/1.-avg_rho)/avg_rho # delta_r=(rho-average rho)/average rho

delta_k_ini=np.array(np.fft.fftn(delta_r_ini))*(L/N)**3. #FT(r to k)=python FT(r to k) * (L/N)^3  # changed from rfftn to fftn


filename='/home/xc298/palmer_scratch/scratch60/reconstruction/output/'+Quijote_ini_dirname+'/large/'+'tsc_mass_ini.dat'
os.makedirs(os.path.dirname(filename), exist_ok=True)
np.save(filename,mass_ini)
np.save('/home/xc298/palmer_scratch/scratch60/reconstruction/output/'\
        +Quijote_ini_dirname+'/large/'+'delta_r_ini.dat',delta_r_ini)
#np.save('/home/xc298/scratch60/reconstruction/output/'\
#        +Quijote_ini_dirname+'/large/'+'delta_k_ini.dat',delta_k_ini)


start=timeit.default_timer() # time it
print("Starting power spectrum for particle_ini")
#k_avg,r_avg,mono_pk_ini,mono_cf_ini,k_counts,r_counts=ps.power_spectrum(N,L,delta_k_ini,pkxibins_spec)
savearr_pk,savearr_xi=ps.power_spectrum(N,L,delta_k_ini,pkxibins_spec)
print("mean delta_r ini=",np.mean(delta_r_ini))
print("Ending power spectrum for particle_ini")


#k_pk=np.array([k_avg,mono_pk_ini,k_counts]).T
#r_cf=np.array([r_avg,mono_cf_ini,r_counts]).T
header2_ini='''
# Simulation path : %s
# Grid size : %i
'''%(sim_path,N)

filename='/home/xc298/project/reconstruction_project/output/'+Quijote_ini_dirname+'/mono_pk_ini.txt'
os.makedirs(os.path.dirname(filename), exist_ok=True)
np.savetxt(filename, savearr_pk, fmt=['%.12f','%.12f','%.12f','%.12f'],header="Monopole power spectrum initial"+header2_ini+\
        '%12s\t%12s\t%12s\t%12s'%('k_avg','P(k) w/ window','P(k) w/o window','k counts'),delimiter='\t')
np.savetxt('/home/xc298/project/reconstruction_project/output/'\
        +Quijote_ini_dirname+'/mono_cf_ini.txt', savearr_xi, fmt=['%.12f','%.12f','%.12f'],header="Monopole correlation function initial"+header2_ini+\
        '%12s\t%12s\t%12s'%('r_avg','xi(r)','r counts'),delimiter='\t')

stop = timeit.default_timer()
print(stop-start)
