from __future__ import print_function
import numpy as np
import networkx as nx

class SIR_node():

    def __init__(self):
        self.set_susceptible()

    def is_infected(self):
        return self.infected

    def is_recovered(self):
        return self.recovered

    def is_susceptible(self):
        return not self.infected and not self.recovered

    def set_infected(self):
        self.infected = True
        self.recovered = False

    def set_recovered(self):
        self.infected = False
        self.recovered = True

    def set_susceptible(self):
        self.infected = False
        self.recovered = False



class Gillespie_SIR():

    def __init__(self,
                 G,              # (dynamic) network. beware! this is going to change during simulation
                 infection_rate, # per S-I link
                 recovery_rate,  # per individual
                 rewiring_rate = 0.,  # per individual
                 infection_seeds = 5,
                 vaccinated = 0,
                 rewire_function = None,
                 mean_degree_function = None
                ):


        if nx.is_directed(G):
            raise TypeError("G is directed but needs to be undirected")

        self.rewire_function = rewire_function
        self.mean_degree_function = mean_degree_function

        if self.rewire_function is None:
            self.rewiring_rate = 0.
        else:
            self.rewiring_rate = rewiring_rate

        self.infection_rate = infection_rate
        self.recovery_rate = recovery_rate

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

            vaccinated_nodes = np.random.choice(list(to_choose_from),vaccinated,replace=False)
            self.recovered.update(vaccinated_nodes)
            #for n in vaccinated_nodes:
            #   self.SIR_nodes[n].set_recovered()


        # add infected seed nodes
        if hasattr(infection_seeds,"__len__"):
            seed_nodes = infection_seeds
        else:
            seed_nodes = np.random.choice(list(self.get_susceptibles()),infection_seeds,replace=False)

        self.infected.update(seed_nodes)
        #for n in vaccinated_nodes:
        #    self.SIR_nodes[n].set_recovered()

        # process new edges for SI links
        self.SI_links = set()
        for newly_inf in self.infected:
            new_edges = [ (newly_inf, n) for n in G.neighbors(newly_inf) ]
            new_SI_links = self.get_SI_links_from_edge_list(new_edges)
            self.SI_links.update(new_SI_links)


    def get_event_rates(self):
        return np.array([ 
                          self.number_of_SI_links() * self.infection_rate,
                          self.number_of_infected() * self.recovery_rate,
                          self.rewiring_rate()      * self.G.number_of_nodes()
                        ])

    def rewire_event(self):

        # get deleted and new edges from rewiring 
        deleted_edges, new_edges = self.rewire_function()
        adjoint_deleted_edges = [ (e[1],e[0]) for e in deleted_edges ]

        # remove all edges from SI links that have been removed
        # due to rewiring
        self.SI_links.remove(deleted_edges)
        self.SI_links.remove(adjoint_deleted_edges)

        # process new edges for SI links
        new_SI_links = self.get_SI_links_from_edge_list(new_edges)
        self.SI_links.update(new_SI_links)

    def get_SI_links_from_edge_list(self,edgelist):
        return [ 
                    e for e in edgelist \
                    if (e[0] in self.infected and e[1] not in self.infected and e[1] not in self.recovered) or \
                       (e[1] in self.infected and e[0] not in self.infected and e[0] not in self.recovered)
               ]

    def recover_event(self):

        recovered = np.random.choice(list(self.infected))
        self.infected.remove(recovered)
        self.recovered.add(recovered)

        deleted_edges = []
        [ deleted_edges.extend([(recovered,n), (n,recovered)]) for n in G.neighbors(recovered) ]
        self.SI_links.remove(deleted_edges)

    def infect_event(self):

        infective_link = np.random.choice(list(self.SI_links))

        if infective_link[0] not in self.infected:
            newly_inf = infective_link[0]
            self.infected.add(infective_link[0])            
        elif infective_link[1] not in self.infected:
            newly_inf = infective_link[1]
            self.infected.add(infective_link[1])
        else:
            raise ValueError("There was a non SI-link in the array of SI links")

        # process new edges for SI links
        new_edges = [ (newly_inf, n) for n in G.neighbors[newly_inf] ]
        new_SI_links = self.get_SI_links_from_edge_list(new_edges)
        self.SI_links.update(new_SI_links)

    def choose_event(self):
        pass

    def simulate(self):
        pass

    def get_number_of_SI_links():
        return len(self.SI_links)

    def get_susceptibles(self):
        return (self.nodes - self.infected) - self.recovered

    def number_of_susceptibles(self):
        return G.number_of_nodes() - self.number_of_infected() - self.number_of_recovered()

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
        return self.S / float(self.G.number_of_nodes())

    def i(self):
        return self.I / float(self.G.number_of_nodes())

    def r(self):
        return self.R / float(self.G.number_of_nodes())

            
if __name__=="__main__":
    from flockworks import flockwork

    F = flockwork(0.5,N=100)

    F.equilibrate()

    infection_rate = 0.1
    recovery_rate = 0.01
    rewiring_rate = 0.01


    sim = Gillespie_SIR(F.G,infection_rate,recovery_rate,rewiring_rate,rewire_function=F.rewire)


