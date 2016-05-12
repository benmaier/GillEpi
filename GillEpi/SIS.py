from __future__ import print_function
import random

import numpy as np
import networkx as nx

import GillEpi.SIR as SIR


def cumstd(arr):
    return np.array([ abs(np.std(arr[:i-1])-np.std(arr[:i]))/np.std(arr[:i]) for i in np.arange(3,len(arr)+1) ])

def cumrelerr(arr):
    return np.array([ np.std(arr[:i]) / np.mean(arr[:i]) / np.sqrt(i) for i in np.arange(2,len(arr)+1) ])

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
                 save_everything = False,
                 without_extinction = False,
                ):

        # SIS behaves differently than SIR. To account for this, we introduce another option for
        # equilibration, to start in a state which should naturally be near equilibrium.
        # BEWARE! this also implies that the given Graph G is in equilibrium
        if infection_seeds == "near_equilibrium":

            self.start_near_equilibrium = True

            if mean_degree_function is not None:
                k0 = mean_degree_function(G)
            else:
                k0 = np.mean( G.degree().values() )

            k0_sqrd = np.mean( np.array(G.degree().values())**2 )

            # basic reproduction number
            R0 = infection_rate / recovery_rate * k0

            if R0 > 1.:
                infection_seeds = int(G.number_of_nodes() * (1. - 1./R0))
            else:
                infection_seeds = 5 

            self.time_scale = 1. / float(infection_rate) * k0 / (k0_sqrd - k0)

        else:
            self.start_near_equilibrium = False
            

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

        self.without_extinction = without_extinction

    def _recover_event(self):

        if self.verbose:
            print("============ recover event")

        if not self.without_extinction or \
           (self.without_extinction and self.number_of_infected()>1):
            recovered = random.sample(self.infected,1)[0]
            self.infected.remove(recovered)
            self.SIR_nodes[recovered].set_susceptible()

            deleted_edges = []
            [ deleted_edges.extend([(recovered,n), (n,recovered)]) for n in self.G.neighbors(recovered) ]

            if self.verbose:
                print("deleted",deleted_edges)

            self.SI_links.difference_update(deleted_edges)

    def _event(self):

        tau,event = self._choose_tau_and_event()
        self.t += tau

        if event==0:
            self._infection_event()
            self.s_of_t.append([ self.t, self.s() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==1:
            self._recover_event()
            self.s_of_t.append([ self.t, self.s() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==2:
            self._rewire_event()
            if self.mean_degree is not None:
                self.k_of_t.append([ self.t, self.mean_degree(self.G) ])

    def equilibrate(self):
        if self.start_near_equilibrium:
            # look for the maximum of rates            
            """
            rate = [ self.infection_rate, self.recovery_rate, self.rewiring_rate ]
            min_rate_event = np.argmin(rates)
            min_rate = self.min_event_rates()[min_rate_event]
            if event > 0:
                event_size = self.G.number_of_nodes()
            else:
                event_size = self.number_of_SI_links()
                equilibration_time = 2 * event_size / max_rate 
            """
            rate = self.infection_rate
            #time_scale = self.G.number_of_nodes() / self.infection_rate
            time_scale = 20 * self.time_scale
            eq_time = time_scale
            self.simulate(eq_time)
        else:
            #raise ValueError("equilibration from low numbers of infected not implemented yet.")
            time_scale = 20. / self.infection_rate
            eq_time = time_scale
            self.simulate(eq_time)

    def equilibrate_variance(self,tol):

        while len(self.s_of_t)<10:
            self._event()

        old_std_sqr = np.std(np.array(self.i_of_t)[:,1])
        old_mean = np.mean(np.array(self.i_of_t)[:,1])
        old_len = len(self.i_of_t)
        n = old_len

        print("first eq done")

        while True:
            self._event()
            new_len = len(self.i_of_t)
            #if old_len<new_len:
            if True:
                new_i = self.i_of_t[-1][1]
                new_mean = n * old_mean / (n+1) + new_i
                new_std_sqr = (n-1) * old_std_sqr / n + (new_i - new_mean) * (new_i - old_mean) / n
                diff = abs(old_std_sqr-new_std_sqr)/ new_std_sqr
                print(diff)
                if diff<tol:
                    break
                old_std_sqr = new_std_sqr
                old_mean = new_mean 
            n = n+1

        #while True:
        #    self._event()
        #    new_len = len(self.i_of_t)
        #    #if old_len<new_len:
        #    if True:
        #        new_i = self.i_of_t[-1][1]
        #        new_mean = old_len * old_mean / new_len + new_i
        #        new_std_sqr = (old_len-1) * old_std_sqr / old_len + (new_i - new_mean) * (new_i - old_mean) / old_len
        #        diff = abs(old_std_sqr-new_std_sqr)/ new_std_sq
        #        print(diff)
        #        if diff<tol:
        #            break
        #        old_std_sqr = new_std_sqr
        #        old_mean = new_mean 

    def simulate(self,tmax):

        new_tmax = self.t + tmax
        while self.t < new_tmax and self.number_of_infected()>0.:
            self._event()


            
if __name__=="__main__":
    from flockworks import flockwork
    import pylab as pl
    import seaborn as sns
    import time
    import scipy as sp


    show_eq = False


    start = time.time()
    print("equilibrating...")
    F = flockwork(0.95,N=500,equilibrate_fast=True)
    #if show_eq:
    #    ts = F.equilibrate(deterministic_equilibration=True,get_time_series=True)
    #else:
    #    F.equilibrate()
    end = time.time()
    print("equilibration done, took", end-start,"seconds")


    R0 = 10.0
    recovery_rate = 10.0
    infection_rate = R0 * recovery_rate / F.mean_degree() 
    rewiring_rate = 1.0


    sim = SIS(F.G,
                        infection_rate,
                        recovery_rate,
                        rewiring_rate,
                        infection_seeds = 0.5,
                        rewire_function = F.rewire,
                        mean_degree_function = F.mean_degree,
                        #verbose = True,
                        without_extinction = True
                        )

    start = time.time()
    print("equilibrate SIS")
    #sim.equilibrate_variance(1e-5)
    sim.simulate(100)
    t_eq = sim.t
    end = time.time()
    print("equilibration done, took", end-start,"seconds")
    sim.simulate(tmax=20/rewiring_rate)

    fig,ax = pl.subplots(3,1,figsize=(9,9))

    s = np.array(sim.s_of_t)
    i = np.array(sim.i_of_t)
    r = np.array(sim.r_of_t)
    ax[0].step(s[:,0],s[:,1])
    ax[0].step(r[:,0],r[:,1])
    ax[0].step(i[:,0],i[:,1])
    observable = sp.integrate.cumtrapz(i[:,1])/i[1:,0]
    ax[2].step(i[1:,0],observable)
    #ax[2].step(i[3:,0],cumstd(sp.integrate.cumtrapz(i[:,1])/i[1:,0]) / np.sqrt(np.arange(3,len(i[:,1])) ))
    ax[2].step(i[2:,0],cumrelerr(observable/i[1:,0]))
    ax[2].set_yscale('log')
    ax[0].plot([t_eq,t_eq],[0,1])

    R0, t = sim.get_R0_of_t()

    if show_eq:
        t_eq = ts['t']-ts['t'][-1]
        R0_eq = ts['Mean Degree'] * infection_rate / recovery_rate
        ax[1].step(np.concatenate((t_eq,t)),np.concatenate((R0_eq,R0)))
    else:
        ax[1].step(t,R0)

    pl.show()
