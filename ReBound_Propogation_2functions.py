
#%%
import rebound
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import progressbar

# sim is simulation class from rebound
# filename is the rebound simulation environment setup file
# filename,MJD_model,SMA,ECC,MJD_model,MA_config in radians
# record trajectory, and interested parameters like Sun distance and Jupiter distance
def trajectory_from_rebound(filename,MJD_model,SMA,ECC,INC,RAAN,AOP,MA,**kwargs):
    sim = rebound.Simulation(filename)
    # Name, a, e, i, long. node, arg. peric., mean anomaly
    # Test Particle
    sim.add(primary=sim.particles[0], a=SMA, e=ECC, inc=INC, Omega=RAAN, omega=AOP, M=MA)
    #   plotting orbit
    #    for orbit in sim.calculate_orbits(primary=sim.particles[0]):
    #        print(orbit)
    #        print("M =", orbit.M)
    
    # integration time parameters
    Noutputs, tspan = kwargs['No_output'], kwargs['tspan']
    Long_term_year = tspan # 1 year is 1 year on Earth
    year = 2.*np.pi # One year in units where G=1
    times = np.linspace(0.,Long_term_year*year, Noutputs)

    MJD,a,e,inc,Omega,omega,M = np.zeros((7,Noutputs))
    second_body_dist = np.zeros(Noutputs) # distance to second body
    
    # integrator parameters
    sim.integrator = kwargs['simintegrator']
    sim.dt = kwargs['simdt'] # we're working in units where 1 year = 2*pi
    sim.ri_ias15.min_dt = kwargs['sim_ias15dt'] # ensure that close encounters do not stall the integration
    sim.ri_ias15.epsilon = kwargs['epsilon'] # control accuracy for ias15
    
    sim.move_to_com()        # We always move to the center of momentum frame before an integration
    ps = sim.particles       # ps is now an array of pointers and will change as the simulation runs
    
    #sim.status() # show status in Cartesian coordinate
    # bar = progressbar.ProgressBar(max_value=max(times/year*365.25)) # progress bar
    for i,time in enumerate(times):
        # bar.update(time/year*365.25) # updating bar
        sim.integrate(time)

        # get orbit element from numerical integration in radians wrt first primary body
        MJD[i] = MJD_model+time/year*365.25
        orbit = ps[-1].calculate_orbit(primary=sim.particles[0])
        a[i],e[i],inc[i] = orbit.a, orbit.e, orbit.inc
        Omega[i],omega[i],M[i] = orbit.Omega, orbit.omega, orbit.M
    
        # get TP's distance wrt specified primary body
        orbit_second_body = ps[-1].calculate_orbit(primary=sim.particles[1]) # wrt second body
        second_body_dist[i] = orbit_second_body.d
    
    #    stacking orbital elements in form of [MJD, a, e, inc, Omega, omega, M] in degree
    oe = np.column_stack((MJD,a,e,np.degrees(inc),np.degrees(Omega),np.degrees(omega),np.degrees(M))) # trajectory
    oe_extend = np.column_stack((oe,second_body_dist)) # trajectory and extended targets
    oe_extend_pd=pd.DataFrame(data=oe_extend[0:,0:],index=None,
                              columns=['MJD','SMA(AU)','ECC','INC(deg)','RAAN(deg)','AOP(deg)','MA(deg)','Second_Dist(AU)'])
    return oe_extend_pd


# sim is simulation class from rebound
# filename is the rebound simulation environment setup file
# x,y,z,vx,vy,vz in AU and AU/TU and TU=yr2pi
# record heliocentric trajectory x,y,z,vx,vy,vz
def trajectory_rv_from_rebound(filename,MJD_model,SMA,ECC,INC,RAAN,AOP,MA,**kwargs):
    sim = rebound.Simulation(filename)
    # Name, a, e, i, long. node, arg. peric., mean anomaly
    # Test Particle
    sim.add(primary=sim.particles[0], a=SMA, e=ECC, inc=INC, Omega=RAAN, omega=AOP, M=MA)    
    #   plotting orbit
    #    for orbit in sim.calculate_orbits(primary=sim.particles[0]):
    #        print(orbit)
    #        print("M =", orbit.M)
    
    # integration time parameters
    Noutputs, tspan = kwargs['No_output'], kwargs['tspan']
    Long_term_year = tspan # 1 year is 1 year on Earth
    year = 2.*np.pi # One year in units where G=1
    times = np.linspace(0.,Long_term_year*year, Noutputs)

    MJD,x,y,z,vx,vy,vz = np.zeros((7,Noutputs))
    
    # integrator parameters
    sim.integrator = kwargs['simintegrator']
    sim.dt = kwargs['simdt'] # we're working in units where 1 year = 2*pi
    sim.ri_ias15.min_dt = kwargs['sim_ias15dt'] # ensure that close encounters do not stall the integration
    sim.ri_ias15.epsilon = kwargs['epsilon'] # control accuracy for ias15

    sim.move_to_com()        # We always move to the center of momentum frame before an integration
    ps = sim.particles       # ps is now an array of pointers and will change as the simulation runs
    
    #sim.status() # show status in Cartesian coordinate
    # bar = progressbar.ProgressBar(max_value=max(times/year*365.25)) # progress bar
    for i,time in enumerate(times):
        # bar.update(time/year*365.25) # updating bar
        sim.integrate(time)
        # get orbit element from numerical integration in radians with Sun
        MJD[i] = MJD_model+time/year*365.25
        # cartesian parameters
        ptc = ps[-1]
        x[i],y[i],z[i] = ptc.x, ptc.y, ptc.z
        vx[i],vy[i],vz[i] = ptc.vx, ptc.vy, ptc.vz

    #    stacking orbital elements in form of [MJD, x, y, z, vx, vy, vz] in AU, yr2pi
    rv = np.column_stack((MJD,x,y,z,vx,vy,vz)) # trajectory x y z vx vy vz
    rv_pd=pd.DataFrame(data=rv[0:,0:], index=None, columns=['MJD','x','y','z','vx','vy','vz'])
    return rv_pd



