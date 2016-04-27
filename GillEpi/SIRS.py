from __future__ import print_function
import random

import numpy as np
import networkx as nx

from SIR_node import SIR_node

class SIRS():

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


        if nx.is_directed(G):
            raise TypeError("G is directed but needs to be undirected")

        self.G = G

        self.rewire_function = rewire_function
        self.mean_degree = mean_degree_function

        self.verbose = verbose

        if self.rewire_function is None:
            self.rewiring_rate = 0.
        else:
            self.rewiring_rate = rewiring_rate

        self.infection_rate = infection_rate
        self.recovery_rate = recovery_rate
        self.susceptible_rate = susceptible_rate

        self.rates = np.array( [ infection_rate, recovery_rate, rewiring_rate ] )

        self.infected = set()
        self.recovered = set()
        self.nodes = set(G.nodes())

        self.SIR_nodes = { n:SIR_node() for n in self.nodes }

        # add vaccinated nodes (are in recovered class right from the beginning)
        if vaccinated>0:
            if hasattr(infection_seeds,"__len__"):
                to_choose_from = self.nodes() - set(infection_seeds)
            else:
                to_choose_from = self.nodes()

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
            removed_SI_links, new_SI_links = self.get_removed_and_new_SI_links_from_edge_list(new_edges)
            self.SI_links.update(new_SI_links)
            self.SI_links.difference_update(removed_SI_links)

        self.t = 0.
        self.s_of_t = [ [ 0., self.s() ] ]
        self.i_of_t = [ [ 0., self.i() ] ]
        self.r_of_t = [ [ 0., self.r() ] ]

        if self.mean_degree is not None:
            self.k_of_t = [ [ 0., self.mean_degree(G) ] ]


    def get_event_rates(self):
        return np.array([ 
                          self.number_of_SI_links()  * self.infection_rate,
                          self.number_of_infected()  * self.recovery_rate,
                          self.number_of_recovered() * self.susceptible_rate,
                          self.G.number_of_nodes()   * self.rewiring_rate
                        ],dtype=float)

    """
    def is_susceptible(self,node):
        return (node not in self.infected) and (node not in self.recovered)

    def get_removed_and_new_SI_links_from_edge_list(self,edgelist):
        new_SI = []
        removed_SI = []
        [ 
              new_SI.append(e) \
              if ( (e[0] in self.infected and self.is_susceptible(e[1]) ) or \
                   (e[1] in self.infected and self.is_susceptible(e[0]) ) )  \
              else removed_SI.extend([e,(e[1],e[0])]) \
              for e in edgelist
        ]
        return removed_SI, new_SI
    """

    def get_removed_and_new_SI_links_from_edge_list(self,edgelist):
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

    def recover_event(self):

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

    def infection_event(self):
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
        removed_SI_links,new_SI_links = self.get_removed_and_new_SI_links_from_edge_list(new_edges)

        if self.verbose:
            print("now infected:", self.infected)
            print("potential new_SI_links:", new_edges)
            print("new_SI_links:", new_SI_links)
            print("removed_SI_links:", removed_SI_links)

        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)

    def susceptible_event(self):

        if self.verbose:
            print("============ susceptible event")
            print("recovered: ", self.recovered)
        newly_susceptible = random.sample(self.recovered,1)[0]
        self.recovered.remove(newly_susceptible)
        self.SIR_nodes[newly_susceptible].set_susceptible()

        new_edges = [ (newly_susceptible, n) for n in self.G.neighbors(newly_susceptible) ]
        removed_SI_links, new_SI_links = self.get_removed_and_new_SI_links_from_edge_list(new_edges)

        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)

    def rewire_event(self):

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
        removed_SI_links, new_SI_links = self.get_removed_and_new_SI_links_from_edge_list(new_edges)
        self.SI_links.update(new_SI_links)
        self.SI_links.difference_update(removed_SI_links)

    def choose_tau_and_event(self):
        rates = self.get_event_rates()
        total_rate = rates.sum()
        tau = np.random.exponential(1./total_rate)
        try:
            event = np.random.choice(len(rates),p=rates/total_rate)
        except ValueError as e:
            print("rates:", rates)
            print("total_rate:", total_rate)
            raise ValueError(e)

        return tau, event

    def event(self):

        tau,event = self.choose_tau_and_event()
        self.t += tau

        if event==0:
            self.infection_event()
            self.s_of_t.append([ self.t, self.s() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==1:
            self.recover_event()
            self.r_of_t.append([ self.t, self.r() ])
            self.i_of_t.append([ self.t, self.i() ])
        elif event==2:
            self.susceptible_event()
            self.r_of_t.append([ self.t, self.r() ])
            self.s_of_t.append([ self.t, self.s() ])
        elif event==3:
            self.rewire_event()
            if self.mean_degree is not None:
                self.k_of_t.append([ self.t, self.mean_degree(self.G) ])


    def simulate(self,tmax):

        while self.t <= tmax and self.number_of_infected()>0 and self.number_of_susceptibles()>0:
            self.event()

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

    def S(self):
        return self.number_of_susceptibles()

    def I(self):
        return self.number_of_infected()

    def R(self):
        return self.number_of_recovered()

    def s(self):
        return self.S() / float(self.G.number_of_nodes())

    def i(self):
        return self.I() / float(self.G.number_of_nodes())

    def r(self):
        return self.R() / float(self.G.number_of_nodes())

    def get_t_R0(self):        
        if self.mean_degree is not None:
            k_of_t = np.array(self.k_of_t)
            t, k = k_of_t[:,0], k_of_t[:,1]
            return t, k * self.infection_rate / self.recovery_rate


            
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

    t,R0 = sim.get_t_R0()

    if show_eq:
        t_eq = ts['t']-ts['t'][-1]
        R0_eq = ts['Mean Degree'] * infection_rate / recovery_rate
        ax[1].step(np.concatenate((t_eq,t)),np.concatenate((R0_eq,R0)))
    else:
        ax[1].step(t,R0)

    pl.show()
