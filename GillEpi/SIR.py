from __future__ import print_function
import random
import sys

import numpy as np
import networkx as nx

from GillEpi import SIR_node

class SIR():

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


        if nx.is_directed(G):
            raise TypeError("G is directed but needs to be undirected")

        self.G = G

        self.rewire_function = rewire_function
        self.mean_degree = mean_degree_function

        self.verbose = verbose

        if self.rewire_function is None:
            self.rewiring_rate = 0.
        else:
            self.rewiring_rate = float(rewiring_rate)

        self.infection_rate = float(infection_rate)
        self.recovery_rate = float(recovery_rate)

        self.rates = np.array( [ infection_rate, recovery_rate, rewiring_rate ] )

        self.infected = set()
        self.recovered = set()
        self.nodes = set(G.nodes())

        self.SIR_nodes = { n:SIR_node() for n in self.nodes }

        if not hasattr(infection_seeds,"__len__") and infection_seeds<1:
            infection_seeds = int(self.G.number_of_nodes()*infection_seeds)

        # add vaccinated nodes (are in recovered class right from the beginning)
        if vaccinated>0:
            if hasattr(infection_seeds,"__len__"):
                to_choose_from = self.nodes - set(infection_seeds)
            else:
                to_choose_from = self.nodes

            vaccinated_nodes = random.sample(to_choose_from,vaccinated)
            self.recovered.update(vaccinated_nodes)
            
            for n in vaccinated_nodes:
               self.SIR_nodes[n].set_recovered()


        # add infected seed nodes
        if hasattr(infection_seeds,"__len__"):
            seed_nodes = infection_seeds
            infection_seeds = len(infection_seeds)
        else:
            seed_nodes = random.sample(self.get_susceptibles(),infection_seeds)

        self.infected.update(seed_nodes)
        for n in seed_nodes:
            self.SIR_nodes[n].set_infected()


        # process new edges for SI links
        self.SI_links = set()
        for newly_inf in self.infected:
            new_edges = [ (newly_inf, n) for n in G.neighbors(newly_inf) ]
            removed_SI_links, new_SI_links = self._get_removed_and_new_SI_links_from_edge_list(new_edges)
            self.SI_links.update(new_SI_links)
            self.SI_links.difference_update(removed_SI_links)

        """
        for e in self.SI_links:
            s = ''
            if self.SIR_nodes[e[0]].is_susceptible():
                s += 'S '
            elif self.SIR_nodes[e[0]].is_infected():
                s += 'I '
            elif self.SIR_nodes[e[0]].is_recovered():
                s += 'R '
            else:
                raise ValueError("hwaht")
            if self.SIR_nodes[e[1]].is_susceptible():
                s += 'S '
            elif self.SIR_nodes[e[1]].is_infected():
                s += 'I '
            elif self.SIR_nodes[e[1]].is_recovered():
                s += 'R '
            else:
                raise ValueError("hwaht")
            print(s)
         """
                

        self.t = 0.
        self.t_max = 0.
        self.s_of_t = [ [ 0., self.s() ] ]
        self.i_of_t = [ [ 0., self.i() ] ]
        self.r_of_t = [ [ 0., self.r() ] ]

        if self.mean_degree is not None:
            self.k_of_t = [ [ 0., self.mean_degree(G) ] ]


    def _get_event_rates(self):
        return np.array([ 
                          self.number_of_SI_links() * self.infection_rate,
                          self.number_of_infected() * self.recovery_rate,
                          self.G.number_of_nodes()  * self.rewiring_rate
                        ],dtype=float)

    def _get_removed_and_new_SI_links_from_edge_list(self,edgelist):
        new_SI = []
        removed_SI = []
        [ 
              new_SI.append(e) \
              if ( (self.SIR_nodes[e[0]].is_infected() and self.SIR_nodes[e[1]].is_susceptible() ) or \
                   (self.SIR_nodes[e[1]].is_infected() and self.SIR_nodes[e[0]].is_susceptible() ) ) \
              else removed_SI.extend([e,(e[1],e[0])]) \
              for e in edgelist
        ]
        return removed_SI, new_SI

    def _recover_event(self):

        if self.verbose:
            print("============ recover event")
        recovered = random.sample(self.infected,1)[0]
        self.infected.remove(recovered)
        self.recovered.add(recovered)
        self.SIR_nodes[recovered].set_recovered()

        deleted_edges = []
        [ deleted_edges.extend([(recovered,n), (n,recovered)]) for n in self.G.neighbors(recovered) ]

        if self.verbose:
            print("deleted",deleted_edges)

        self.SI_links.difference_update(deleted_edges)

    def _infection_event(self):
        #if self.verbose:
        if self.verbose:   
            print("============= infection event")
            print("infected:", self.infected)
            print("recovered:", self.recovered)
            print("SI link", self.SI_links)

        infective_link = random.sample(self.SI_links,1)[0]

        if self.verbose:
            print("infective_link:", infective_link)

        if self.SIR_nodes[ infective_link[0] ].is_infected() and self.SIR_nodes[ infective_link[1] ].is_susceptible():
            newly_inf = infective_link[1]
            self.infected.add(infective_link[1])
        elif self.SIR_nodes[ infective_link[1] ].is_infected() and self.SIR_nodes[ infective_link[0] ].is_susceptible():
            newly_inf = infective_link[0]
            self.infected.add(infective_link[0])
        else:
            raise ValueError("There was a non SI-link in the array of SI links. This shouldn't happen.")        

        self.SIR_nodes[newly_inf].set_infected()

        # process new edges for SI links
        new_edges = [ (newly_inf, n) for n in self.G.neighbors(newly_inf) ]
        removed_SI_links,new_SI_links = self._get_removed_and_new_SI_links_from_edge_list(new_edges)

        if self.verbose:
            print("now infected:", self.infected)
            print("potential new_SI_links:", new_edges)
            print("new_SI_links:", new_SI_links)
            print("removed_SI_links:", removed_SI_links)

        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)


    def _rewire_event(self):

        # get deleted and new edges from rewiring 
        deleted_edges, new_edges = self.rewire_function()

        if self.verbose:
            print("============= rewiring event")
            print("rewired:", deleted_edges,new_edges)

        adjoint_deleted_edges = [ (e[1],e[0]) for e in deleted_edges ]

        # remove all edges from SI links that have been removed
        # due to rewiring
        self.SI_links.difference_update(deleted_edges)
        self.SI_links.difference_update(adjoint_deleted_edges)

        # process new edges for SI links
        removed_SI_links, new_SI_links = self._get_removed_and_new_SI_links_from_edge_list(new_edges)
        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)

    def _choose_tau_and_event(self):
        rates = self._get_event_rates()
        #print(self.G.number_of_edges(),self.number_of_SI_links())
        #print(self.number_of_SI_links(),self.number_of_infected()*(self.G.number_of_nodes()-self.number_of_infected()))
        total_rate = rates.sum()
        tau = np.random.exponential(1./total_rate)
        try:
            event = np.random.choice(len(rates),p=rates/total_rate)
        except ValueError as e:
            print("rates:", rates)
            print("total_rate:", total_rate)
            raise ValueError(e)

        return tau, event

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
            self._rewire_event()
            if self.mean_degree is not None:
                self.k_of_t.append([ self.t, self.mean_degree(self.G) ])


    def simulate(self):

        while self.number_of_infected() > 0 and self.number_of_susceptibles() > 0:
            self._event()

    def number_of_SI_links(self):
        return len(self.SI_links)

    def get_susceptibles(self):
        return (self.nodes - self.infected) - self.recovered

    def number_of_susceptibles(self):
        return self.G.number_of_nodes() - self.number_of_infected() - self.number_of_recovered()

    def number_of_infected(self):
        return len(self.infected)

    def number_of_recovered(self):
        return len(self.recovered)

    def s(self,normed=True):
        if normed:
            return self.number_of_susceptibles() / float(self.G.number_of_nodes())
        else:
            return self.number_of_susceptibles()

    def i(self,normed=True):
        if normed:
            return self.number_of_infected() / float(self.G.number_of_nodes())
        else:
            return self.number_of_infected()

    def r(self,normed=True):
        if normed:
            return self.number_of_recovered() / float(self.G.number_of_nodes())
        else:
            return self.number_of_recovered()

    def get_outbreak_size(self,normed=True):
        """return the size of the outbreak (number of infected and recovered)"""
        if normed:
            return 1. - self.s()
        else:
            return self.G.number_of_nodes() - self.number_of_susceptibles()

    def _get_max_t(self):
        """return the time of the last event"""
        """
        if hasattr(self,'k_of_t'):
            return max([ 
                            self.s_of_t[-1][0],
                            self.i_of_t[-1][0],
                            self.r_of_t[-1][0],
                            self.k_of_t[-1][0],
                      ])
        else:
            return max([ 
                            self.s_of_t[-1][0],
                            self.i_of_t[-1][0],
                            self.r_of_t[-1][0],
                      ])
        """
        return self.t_max

    def _get_x_of_t(self,arr,normed=True):
        """get the time of the last event, append it to the list and pass back an nd.array"""
        t_max = self._get_max_t()
        arr = list(arr)

        if arr[-1][0]<t_max:
            arr.append([t_max,arr[-1][1]])

        arr = np.array(arr)
        if normed:
            return arr[:,1], arr[:,0]
        else:
            return arr[:,1]*self.G.number_of_nodes(), arr[:,0]


    def get_s_of_t(self,normed=True):
        return self._get_x_of_t(self.s_of_t,normed)

    def get_i_of_t(self,normed=True):
        return self._get_x_of_t(self.i_of_t,normed)

    def get_r_of_t(self,normed=True):
        return self._get_x_of_t(self.r_of_t,normed)

    def get_k_of_t(self):
        if self.mean_degree is not None:
            return self._get_x_of_t(self.k_of_t,normed=True)
        else:
            raise ValueError("degree has not been calculated since a function was not provided")

    def get_R0_of_t(self):        
        k,t = self.get_k_of_t()
        return k * self.infection_rate / self.recovery_rate, t


    def _get_x_starting_at_t(self,arr,t_start,normed=True):
        x,t = self._get_x_of_t(arr,normed=True)
        ndcs = np.where(t>=t_start)
        return x[ndcs], t[ndcs]

    def get_s_starting_at_t(self,t_start,normed=True):
        return self._get_x_starting_at_t(self.s_of_t,t_start,normed)

    def get_i_starting_at_t(self,t_start,normed=True):
        return self._get_x_starting_at_t(self.i_of_t,t_start,normed)

    def get_r_starting_at_t(self,t_start,normed=True):
        return self._get_x_starting_at_t(self.r_of_t,t_start,normed)

    def get_k_starting_at_t(self,t_start):
        if self.mean_degree is not None:
            return self._get_x_starting_at_t(self.k_of_t,t_start,normed=True)
        else:
            raise ValueError("degree has not been calculated since a function was not provided")

    def get_R0_starting_at_t(self,t_start):        
        k,t = self.get_k_starting_at_t(t_start)
        return k * self.infection_rate / self.recovery_rate, t


            
if __name__=="__main__":
    from flockworks import flockwork
    import pylab as pl
    import seaborn as sns
    import time


    show_eq = False

    F = flockwork(0.8,N=100)

    start = time.time()
    print("equilibrating...")
    if show_eq:
        ts = F.equilibrate(deterministic_equilibration=True,get_time_series=True)
    else:
        F.equilibrate()
    end = time.time()
    print("equilibration done, took", end-start,"seconds")

    infection_rate = 1.
    recovery_rate = 1.
    rewiring_rate = 1.


    sim = SIR(F.G,
                        infection_rate,
                        recovery_rate,
                        rewiring_rate,
                        infection_seeds = 5,
                        rewire_function = F.rewire,
                        mean_degree_function = F.mean_degree,
                        #verbose = True,
                        )

    sim.simulate()

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
