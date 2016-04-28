from __future__ import print_function
import random

import numpy as np
import networkx as nx

from GillEpi.SIR_node import SIR_node
from GillEpi.SIR import SIR

class SIS(SIR):

    def __init__(self,
                 G,              # (dynamic) network. beware! this is going to change during simulation
                 infection_rate, # per S-I link
                 recovery_rate,  # per individual
                 rewiring_rate = 0.,  # per individual
                 infection_seeds = 5,
                 vaccinated = 0,
                 rewire_function = None,
                 mean_degree_function = None,
                 verbose = False,
                 save_everything = False
                ):

        SIR.__init__(self,
                 G,              # (dynamic) network. beware! this is going to change during simulation
                 infection_rate, # per S-I link
                 recovery_rate,  # per individual
                 rewiring_rate = rewiring_rate,  # per individual
                 infection_seeds = infection_seeds,
                 vaccinated = vaccinated,
                 rewire_function = rewire_function,
                 mean_degree_function = mean_degree_function,
                 verbose = verbose,
                 save_everything = save_everything
                )

    def recover_event(self):

        if self.verbose:
            print("============ recover event")
        recovered = random.sample(self.infected,1)[0]
        self.infected.remove(recovered)
        self.SIR_nodes[recovered].set_susceptible()

        deleted_edges = []
        [ deleted_edges.extend([(recovered,n), (n,recovered)]) for n in self.G.neighbors(recovered) ]

        if self.verbose:
            print("deleted",deleted_edges)

        self.SI_links.difference_update(deleted_edges)

    def event(self):

        tau,event = self.choose_tau_and_event()
        self.t += tau

        if event==0:
            self.infection_event()
            self.s_of_t.append([ self.t, self.s() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==1:
            self.recover_event()
            self.s_of_t.append([ self.t, self.s() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==2:
            self.rewire_event()
            if self.mean_degree is not None:
                self.k_of_t.append([ self.t, self.mean_degree(self.G) ])


    def simulate(self,tmax):

        while self.t < tmax:
            self.event()


            
if __name__=="__main__":
    from flockworks import flockwork
    import pylab as pl
    import seaborn as sns
    import time


    show_eq = False

    F = flockwork(0.7,N=1000)

    start = time.time()
    print("equilibrating...")
    if show_eq:
        ts = F.equilibrate(deterministic_equilibration=True,get_time_series=True)
    else:
        F.equilibrate()
    end = time.time()
    print("equilibration done, took", end-start,"seconds")

    infection_rate = 1.5 
    recovery_rate = 1.
    rewiring_rate = 1.


    sim = SIS(F.G,
                        infection_rate,
                        recovery_rate,
                        rewiring_rate,
                        infection_seeds = 5,
                        rewire_function = F.rewire,
                        mean_degree_function = F.mean_degree,
                        #verbose = True,
                        )

    sim.simulate(tmax=10)

    fig,ax = pl.subplots(2,1,figsize=(9,6))

    s = np.array(sim.s_of_t)
    i = np.array(sim.i_of_t)
    r = np.array(sim.r_of_t)
    ax[0].step(s[:,0],s[:,1])
    ax[0].step(r[:,0],r[:,1])
    ax[0].step(i[:,0],i[:,1])

    t,R0 = sim.get_t_R0()

    if show_eq:
        t_eq = ts['t']-ts['t'][-1]
        R0_eq = ts['Mean Degree'] * infection_rate / recovery_rate
        ax[1].step(np.concatenate((t_eq,t)),np.concatenate((R0_eq,R0)))
    else:
        ax[1].step(t,R0)

    pl.show()
