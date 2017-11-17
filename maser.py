import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# coupled differential equations
#\frac{d\Delta E(t)}{dt} = -q_e F v_e \frac{1}{2} sin(\phi(t))\\
#\frac{d\phi(t)}{dt} = \frac{qB}{m_e^2}(\Delta E(t)) 
# params = (q_e F v_e \frac{1}{2} ,  \frac{qB}{m_e^2} ) 

def f(y,t,params):
    deltaE, phi = y
    Econst, phiconst = params
    return  [Econst*np.sin(phi), -phiconst*deltaE]

def feloss(y,t,params):
    deltaE, phi = y
    Econst, phiconst,Elossconstant = params
    return  [Econst*np.sin(phi) - Elossconstant, -phiconst*deltaE]

figuresdirectory = "writeup/"

# if we we use constants such that q=1, F = V/m, v = m/s then dE/dt = eV/s
# if we use constants such that qB/m = 30 GHz, m = 511000 eV, then dphi/dt is per second
me=510998.0
e0 = 18570.0
ve = 2.99e8*np.sqrt(1-(me/(me+e0))**2)
qe = 1.0
F = 0.1e0
wcyc = 100e9

if False:
    params = [qe*F*ve/2, 160e9/me]

    plt.clf()
    fig = plt.figure(1, figsize=(16,16))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_ylim(-20,20)
    ax1.set_xlim(-np.pi,np.pi)

    ax2.set_ylim(-20,20)


    for i,e0 in enumerate(np.arange(-20.5,20,1.0)):
        initialconditions = [e0, 0.0]
        times = np.arange(0.,1e-4,2e-8)

        psoln = odeint(f,initialconditions,times,args=(params,))
        # force phases to wrap
        phases = psoln[:,1]
        phases = np.mod(phases+np.pi,2*np.pi) - np.pi
        energies = psoln[:,0]
    
        ax1.scatter(phases,energies,s=1)
        print i, i%4
        if (i%3 != 0):
            ax2.scatter(times[0:1]*1e6,energies[0:1],s=1)
        else:
            ax2.scatter(times[0:300]*1e6,energies[0:300],s=1)

        ax1.set_xlabel("phi")
        ax1.set_ylabel("delta E")

        ax2.set_xlabel("time (usec)")
        ax2.set_ylabel("delta E")
        fig.savefig(figuresdirectory+"maser_phase_space_map.png")
        plt.show()

#######
# another version where we've added 1 fW of radiative energy loss
# 1 fW = 1e-15 J/s, 1 J = 1.6e-19 eV
if False:
    F = 0.01
    params = [qe*F*ve/2, 160e9/me, 1e-15/1.6e-19]

    plt.clf()
    fig = plt.figure(1, figsize=(16,16))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_ylim(-20,20)
    ax1.set_xlim(-np.pi,np.pi)

    ax2.set_ylim(-20,20)


    for i,e0 in enumerate(np.arange(-20.,20,1.0)):
        initialconditions = [e0, 0.0]
        times = np.arange(0.,3e-4,2e-8)

        psoln = odeint(feloss,initialconditions,times,args=(params,))
        # force phases to wrap
        phases = psoln[:,1]
        phases = np.mod(phases+np.pi,2*np.pi) - np.pi
        energies = psoln[:,0]
    
        ax1.scatter(phases,energies,s=1)
        print i, i%4, e0
        if (i%3 != 0 and abs(e0) > 2.2 ):
            ax2.scatter(times[0:1]*1e6,energies[0:1],s=1)
        else:
            ax2.scatter(times*1e6,energies,s=1)

      
        ax1.set_xlabel("phi")
        ax1.set_ylabel("delta E")

        ax2.set_xlabel("time (usec)")
        ax2.set_ylabel("delta E")
        fig.savefig(figuresdirectory+"maserEL_phase_space_map.png")
        plt.show()


#######
# another version where we've added 1 fW of radiative energy loss
# 1 fW = 1e-15 J/s, 1 J = 1.6e-19 eV
if False:
    F = 0.01
    params = [qe*F*ve/2, 160e9/me, 1e-15/1.6e-19]

    plt.clf()
    fig = plt.figure(1, figsize=(16,16))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.set_ylim(-5,5)
    ax1.set_xlim(-np.pi,np.pi)

    ax2.set_ylim(-5,5)


    for i,e0 in enumerate([3.5,2.0,1.0,0.25]):
        initialconditions = [e0, 0.0]
        times = np.arange(0.,6e-4,2e-8)

        psoln = odeint(feloss,initialconditions,times,args=(params,))
        # force phases to wrap
        phases = psoln[:,1]
        phases = np.mod(phases+np.pi,2*np.pi) - np.pi
        energies = psoln[:,0]
    
        ax1.scatter(phases,energies,s=1)
        print i, i%4, e0
        ax2.scatter(times*1e6,energies,s=1)


        ax1.set_xlabel("phi")
        ax1.set_ylabel("delta E")

        ax2.set_xlabel("time (usec)")
        ax2.set_ylabel("delta E")
        fig.savefig(figuresdirectory+"maser_dodge_phase_space_map.png")
        plt.show()

        
if True:
    plt.clf()

    fig = plt.figure(1, figsize=(16,16))
    ax3 = plt.subplot(313)
    ax1 = plt.subplot(311,sharex=ax3)
    ax2 = plt.subplot(312,sharex=ax3)
#    ax4 = plt.subplot(224)

    minor_ticks = np.arange(-10, 10, 1)  
    
    for F in (0.001,0.00333,0.01):
        params = [qe*F*ve/2, 160e9/me]


        energylist = np.arange(-10.,10,0.2505)
        erange = np.zeros_like(energylist)
        period = np.zeros_like(energylist)
        powers = np.zeros_like(energylist)
        for i,e0 in enumerate(energylist):
            initialconditions = [e0, 0.0]
            times = np.arange(0.,1e-4,2e-9)

            psoln = odeint(f,initialconditions,times,args=(params,))
            erange[i] = max(psoln[:,0])-min(psoln[:,0])
            phases = psoln[:,1]
            phases = np.mod(phases+np.pi,2*np.pi) - np.pi
            # cute code for a list of zero crossings in the phase array
            # 2nd zero crossing is the first full cycle
            myzeros = np.where(np.diff(np.sign(phases)))[0]
            if (len(myzeros) < 3):
                period[i] = len(phases)*(times[1]-times[0])
            else:
                period[i] = myzeros[2]*(times[1]-times[0])
            
        powers = (2*erange/period)*1.6e-19*1e15
        
        ax1.plot(energylist,erange,label="F=%.0e V/m"%F)
        ax2.semilogy(energylist,1.0/period,label="F=%.0e V/m"%F)
        ax3.semilogy(energylist,powers,label="F=%.0e V/m"%F)
#        ax4.scatter(1.0/period,powers,label="F=%.0e V/m"%F)

        
#    ax4.set_xscale('log')
#    ax4.set_yscale('log')

    ax1.set_xticks(minor_ticks, minor=True)
    ax2.set_xticks(minor_ticks, minor=True)
    ax3.set_xticks(minor_ticks, minor=True)


    ax1.grid(which='minor', alpha=0.2)
    ax2.grid(which='minor', alpha=0.2)
    ax3.grid(which='minor', alpha=0.2)
    
    ax3.set_xlabel("energy (eV)")
    ax1.legend(loc="upper right")
    ax1.set_ylabel("energy range")
    ax2.set_ylabel("orbit frequency")
    ax3.set_ylabel("power (fW)")
    ax3.set_ylim(1.0,100)
    plt.savefig(figuresdirectory+"range_freq_power.png")

    plt.show()
