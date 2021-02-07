
#%%
import rebound
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
from ReBound_Propogation_2functions import trajectory_from_rebound
from ReBound_Propogation_2functions import trajectory_rv_from_rebound
from matplotlib import cm
import multiprocessing
from joblib import Parallel, delayed
import time,os
import datetime, julian
from cislunar_constant import *

folder = "../Data/"
# date
date = "2010-03-26 02:05"
date_info = datetime.datetime.strptime(date,"%Y-%m-%d %H:%M")
MJD_model = julian.to_jd(date_info,fmt='mjd')
# rebound simulation
# sim = rebound.Simulation()
# sim.add('Geocenter',date = date)
# sim.add('Luna')
# sim.add('Sun')

# save to binary setup file
simfile = f"CislunarSun_{MJD_model}.bin"
# sim.save(simfile)
sim = rebound.Simulation(simfile)
print(f'{sim.units}')
start_time = time.time()

earth = sim.particles[0]
moon = sim.particles[1]
moon_orbit = moon.calculate_orbit(primary=sim.particles[0]) # j2000 
sun = sim.particles[2]
sun_orbit = sun.calculate_orbit(primary=sim.particles[0])
print(moon_orbit,sun_orbit)

# create meshgrid from self defined values
SMA_min, SMA_max = 20 * Re / AU, 80 * Re / AU #AU
ECC_min, ECC_max = 0.00001, 0.99999
No_steps_SMA = 100
No_steps_ECC = 20
SMA_linspace = np.linspace(SMA_min,SMA_max,No_steps_SMA) 
ECC_linspace = np.linspace(ECC_min,ECC_max,No_steps_ECC) 
SMA_grid,ECC_grid = np.meshgrid(SMA_linspace,ECC_linspace)

# flatten for calculation
SMA_flatlist = SMA_grid.flatten().tolist()
ECC_flatlist = ECC_grid.flatten().tolist()
print(len(SMA_flatlist))

# take same INC,RAAN,AOP,MA in rad
INC_moon,RAAN_moon,AOP_moon,MA_moon = moon_orbit.inc,moon_orbit.Omega,moon_orbit.omega,moon_orbit.M

# for Trojan case only
MA_config = -60
# for Greek case only
# MA_config = 60
# for custom case
#MA_config = 0

# custom INC,RAAN,AOP,MA
INC,RAAN,AOP= INC_moon,RAAN_moon,AOP_moon
MA = MA_moon + np.deg2rad(MA_config)

# rebound propagation detail
simintegrator = 'ias15' #'mercurius'
simdt = 0.0001 # 2pi is 1 year, time follows equation t = simdt/(2*pi)*1*365.25 (day)
sim_ias15dt = 0 # ensure that close encounters do not stall the integration
epsilon = 1e-8 # control accuracy for ias15
tspan = 9/12 # 1 year is 1 year on Earth 365.25 solar days/earth days
No_output = 1000
rb_detail = {'simintegrator':simintegrator, 'simdt':simdt, 'sim_ias15dt':sim_ias15dt,
            'epsilon':epsilon, 'No_output':No_output, 'tspan':tspan}

# Count number of cores
numCores = multiprocessing.cpu_count()
# numCores = 11

# trajectory write function
def trajectory_write(folder,simfile,MJD,SMA,ECC,INC,RAAN,AOP,MA,rb_detail):
    print(SMA)
    oe_extend_pd = trajectory_from_rebound(simfile,MJD,SMA,ECC,INC,RAAN,AOP,MA,**rb_detail)
    # rv_pd        = trajectory_rv_from_rebound(simfile,MJD,SMA,ECC,INC,RAAN,AOP,MA)
    filename = simfile[0:-4] # delet string '.bin'
    trajectory_foldername = folder + f"Trajectory_TrojanCase_{filename}_{No_steps_SMA}x{No_steps_ECC}_{simintegrator}"
    if not os.path.exists(trajectory_foldername): # create trajectory folder if not existed
        os.makedirs(trajectory_foldername)
    
    trajectory_filename="trajectory_"+"SMA"+str(SMA)+"_"+"ECC"+str(ECC)+"_oe.csv"
    # trajectory_filename="trajectory_"+"SMA"+str(SMA)+"_"+"ECC"+str(ECC)+"_oe.h5"
    oe_extend_pd.to_csv(trajectory_foldername+"/"+trajectory_filename)
    # oe_extend_pd.to_hdf(trajectory_foldername+"/"+trajectory_filename, key = 'trajectory', mode='w', format='table')

    # trajectory_rv_filename="trajectory_"+"SMA"+str(SMA)+"_"+"ECC"+str(ECC)+"_rv.csv"
#    trajectory_rv_filename="trajectory_"+"SMA"+str(SMA)+"_"+"ECC"+str(ECC)+"_rv.h5"
    # rv_pd.to_csv(trajectory_foldername+"/"+trajectory_rv_filename)
#    rv_pd.to_hdf(trajectory_foldername+"/"+trajectory_rv_filename, key = 'trajectory', mode='w', format='table')
    return

# Run processes trajectory write functions
Parallel(n_jobs = numCores)(delayed(trajectory_write)(folder,simfile,MJD_model,SMA,ECC,INC,RAAN,AOP,MA,rb_detail) for SMA,ECC in zip(SMA_flatlist,ECC_flatlist))

# time
print(time.time()-start_time)
