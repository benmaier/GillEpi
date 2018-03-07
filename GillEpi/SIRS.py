from __future__ import print_function
import random

import numpy as np
import networkx as nx

import GillEpi.SIR as SIR

class SIRS(SIR):

    def __init__(self,
                 G,              # (dynamic) network. beware! this is going to change during simulation
                 infection_rate, # per S-I link
                 recovery_rate,  # per individual
                 susceptible_rate,  # per individual
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

        self.susceptible_rate = susceptible_rate

    def _get_event_rates(self):
        return np.array([ 
                          self.number_of_SI_links()  * self.infection_rate,
                          self.number_of_infected()  * self.recovery_rate,
                          self.number_of_recovered() * self.susceptible_rate,
                          self.G.number_of_nodes()   * self.rewiring_rate
                        ],dtype=float)

    def _susceptible_event(self):

        if self.verbose:
            print("============ susceptible event")
            print("recovered: ", self.recovered)
        newly_susceptible = random.sample(self.recovered,1)[0]
        self.recovered.remove(newly_susceptible)
        self.SIR_nodes[newly_susceptible].set_susceptible()

        new_edges = [ (newly_susceptible, n) for n in self.G.neighbors(newly_susceptible) ]
        removed_SI_links, new_SI_links = self._get_removed_and_new_SI_links_from_edge_list(new_edges)

        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)

    def _event(self):

        tau,event = self._choose_tau_and_event()
        self.t += tau
        self.t_max = self.t

        if event==0:
            self._infection_event()
            self.s_of_t.append([ self.t, self.s() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==1:
            self._recover_event()
            self.r_of_t.append([ self.t, self.r() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==2:
            self._susceptible_event()
            self.r_of_t.append([ self.t, self.r() ])
            self.s_of_t.append([ self.t, self.s() ])
        elif event==3:
            self._rewire_event()
            if self.mean_degree is not None:
                self.k_of_t.append([ self.t, self.mean_degree(self.G) ])


    def simulate(self,tmax):

        while self.t <= tmax:
            self._event()


            
if __name__=="__main__":
    from flockworks import flockwork
    import pylab as pl
    import seaborn as sns
    import time


    show_eq = False

    F = flockwork(0.7,N=100)

    start = time.time()
    print("equilibrating...")
    if show_eq:
        ts = F.equilibrate(deterministic_equilibration=True,get_time_series=True)
    else:
        F.equilibrate()
    end = time.time()
    print("equilibration done, took", end-start,"seconds")

    infection_rate = 2. 
    recovery_rate = 1.
    susceptible_rate = 1.
    rewiring_rate = 1.


    sim = SIRS(F.G,
                        infection_rate,
                        recovery_rate,
                        susceptible_rate,
                        rewiring_rate,
                        infection_seeds = 2,
                        rewire_function = F.rewire,
                        mean_degree_function = F.mean_degree,
                        #verbose = True,
                        )

    sim.simulate(40)

    fig,ax = pl.subplots(2,1,figsize=(9,6))

    s = np.array(sim.s_of_t)
    i = np.array(sim.i_of_t)
    r = np.array(sim.r_of_t)
    ax[0].step(s[:,0],s[:,1])
    ax[0].step(r[:,0],r[:,1])
    ax[0].step(i[:,0],i[:,1])

    R0, t = sim.get_R0_of_t()

    if show_eq:
        t_eq = ts['t']-ts['t'][-1]
        R0_eq = ts['Mean Degree'] * infection_rate / recovery_rate
        ax[1].step(np.concatenate((t_eq,t)),np.concatenate((R0_eq,R0)))
    else:
        ax[1].step(t,R0)

    pl.show()
