from numpy import *
import networkx as nx
#import numpy.random as random
#import matplotlib
#matplotlib.use('TkAgg')
#import pylab as pl
#from matplotlib import animation

class SIR:

    def __init__(self,G,infection_probability,\
                        recovery_probability,\
                        outbreak_probability = 0.0,\
                        save_node_changes = False,
                        patients_zero = None):
        
        self.DONE = False

        self.save_node_changes = save_node_changes
        self.r_I = infection_probability
        self.r_R = recovery_probability
        self.r_O = outbreak_probability
        self.G = G
        self.patients_zero = patients_zero

        self.newly_infected = []
        self.newly_susceptible = []
        self.newly_recovered = []
        self.activated_edges = []

        self.is_susceptible = {}
        self.infected = set([])
        self.susceptible = set(G.nodes())

        for n in self.susceptible:
            self.is_susceptible[n] = True

        self.I = [0]
        self.S = [len(G)]
        self.T = 0

        
        if not patients_zero is None:
            self.I[-1] += len(patients_zero)
            self.S[-1] -= len(patients_zero)
            for p in patients_zero:
                self.is_susceptible[p] = False
            self.susceptible = self.susceptible - set(patients_zero)
            self.infected = self.infected | set(patients_zero)
            if save_node_changes:
                self.newly_infected.append(patients_zero)

    def reset(self):
        self.__init__(self.G, self.r_I, self.r_R, self.r_O, self.save_node_changes,self.patients_zero)

    def step(self,N_t=1):

        if self.DONE:
            print("No susceptible nodes. Simulating further will not have any effect")
            return 

        for t in range(N_t):
            
            infected_from_spreading, recovered, active_edges = self.spread_and_recover()
            if self.r_O > 0.:
                infected_from_outbreak = self.outbreak()
            else:
                infected_from_outbreak = set([])

            inf = infected_from_spreading | infected_from_outbreak
            dI = len(inf) - len(recovered)
            dS = (-1) * len(inf)

            if self.save_node_changes:
                self.newly_infected.append(inf) 
                self.newly_susceptible.append(set([]))
                self.newly_recovered.append(recovered)
                self.activated_edges.append(active_edges)

            self.I.append(self.I[-1] + dI)
            self.S.append(self.S[-1] + dS)

            self.T+=1

            if self.I[-1] == 0:
                self.DONE = True
                break
        
    def outbreak(self):
        to_be_infected = set()
        for s in self.susceptible:
            if random.rand() < self.r_O:
                to_be_infected.add(s)

        self.infected = self.infected | to_be_infected
        self.susceptible = self.susceptible - to_be_infected

        for t in to_be_infected:
            self.is_susceptible[t] = False

        return to_be_infected    

    def spread_and_recover(self):

        to_be_infected = set()
        to_be_recovered = set()
        active_edges = set()

        for i in self.infected:
            #for neigh in G.neighbors(i):
            for e in self.G.edges(i):
                neigh = e[1]
                if self.is_susceptible[neigh] and random.rand() < self.r_I:
                    to_be_infected.add(neigh)
                    active_edges.add(e)
                    
            if random.rand() < self.r_R:
                to_be_recovered.add(i)

        self.susceptible = self.susceptible - to_be_infected
        self.infected = self.infected | to_be_infected
        self.infected = self.infected - to_be_recovered

        for t in to_be_infected:
            self.is_susceptible[t] = False
                
        return to_be_infected, to_be_recovered, active_edges


class SIS:

    def __init__(self,G,infection_probability,\
                        recovery_probability,\
                        outbreak_probability = 0.0,\
                        save_node_changes = False,
                        patients_zero = None):
        
        self.DONE = False

        self.save_node_changes = save_node_changes
        self.r_I = infection_probability
        self.r_R = recovery_probability
        self.r_O = outbreak_probability
        self.G = G
        self.patients_zero = patients_zero

        self.newly_infected = []
        self.newly_susceptible = []
        self.newly_recovered = []
        self.activated_edges = []

        self.is_susceptible = {}
        self.infected = set([])
        self.susceptible = set(G.nodes())

        for n in self.susceptible:
            self.is_susceptible[n] = True

        self.I = [0]
        self.S = [len(G)]
        self.T = 0

        
        if not patients_zero is None:
            self.I[-1] += len(patients_zero)
            self.S[-1] -= len(patients_zero)
            for p in patients_zero:
                self.is_susceptible[p] = False
            self.susceptible = self.susceptible - set(patients_zero)
            self.infected = self.infected | set(patients_zero)
            if save_node_changes:
                self.newly_infected.append(patients_zero)

    def reset(self):
        self.__init__(self.G, self.r_I, self.r_R, self.r_O, self.save_node_changes,self.patients_zero)

    def step(self,N_t=1):

        if self.DONE:
            print("No susceptible nodes. Simulating further will not have any effect")
            return 

        for t in range(N_t):
            
            infected_from_spreading, recovered, active_edges = self.spread_and_recover()
            if self.r_O > 0.:
                infected_from_outbreak = self.outbreak()
            else:
                infected_from_outbreak = set([])

            inf = infected_from_spreading | infected_from_outbreak
            dI = len(inf) - len(recovered)
            dS = (-1) * len(inf) + len(recovered)

            if self.save_node_changes:
                self.newly_infected.append(inf) 
                self.newly_susceptible.append(recovered)
                self.newly_recovered.append(set([]))
                self.activated_edges.append(active_edges)

            self.I.append(self.I[-1] + dI)
            self.S.append(self.S[-1] + dS)

            self.T+=1

            if self.I[-1] == 0:
                self.DONE = True
                break
        
    def outbreak(self):
        to_be_infected = set()
        for s in self.susceptible:
            if random.rand() < self.r_O:
                to_be_infected.add(s)

        self.infected = self.infected | to_be_infected
        self.susceptible = self.susceptible - to_be_infected

        for t in to_be_infected:
            self.is_susceptible[t] = False

        return to_be_infected    

    def spread_and_recover(self):

        to_be_infected = set()
        to_be_recovered = set()
        active_edges = set()

        for i in self.infected:
            #for neigh in G.neighbors(i):
            for e in self.G.edges(i):
                neigh = e[1]
                if self.is_susceptible[neigh] and random.rand() < self.r_I:
                    to_be_infected.add(neigh)
                    active_edges.add(e)
                    
            if random.rand() < self.r_R:
                to_be_recovered.add(i)

        self.susceptible = self.susceptible - to_be_infected
        self.susceptible = self.susceptible | to_be_recovered
        self.infected = self.infected | to_be_infected
        self.infected = self.infected - to_be_recovered

        for t in to_be_infected:
            self.is_susceptible[t] = False

        for t in to_be_recovered:
            self.is_susceptible[t] = True
                
        return to_be_infected, to_be_recovered, active_edges

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('TkAgg')
    import pylab as pl
    from matplotlib import animation
    from networkx.drawing.nx_agraph import graphviz_layout

    """
    #MHR parameters
    B = 3
    L = 4
    p = 1.0
    xi = 0.6

    N=float(B**L)
        
    seed = 476 #no other connection
    #simulation parameters
    infection_probability = 0.38
    recovery_probability = 0.15
    outbreak_probability = 1e-5
    N_t = 200
    """
    """
    #MHR parameters
    B = 8
    L = 3
    p = 0.9
    xi = 0.4

    N=float(B**L)
        
    seed = 76 #no other connection
    #simulation parameters
    infection_probability = 0.01
    recovery_probability = 0.001
    outbreak_probability = 0
    N_t = 200
    """
    #MHR parameters
    B = 6
    L = 3
    k = 7
    xi = 0.3
    if xi<1. or xi>1.:
        p = k/float(B-1)*(1-xi)/(1-xi**L)
    else:    
        p = k/float(B-1)/L

    N=float(B**L)
        
    seed = 46 #no other connection
    #simulation parameters
    infection_probability = 0.2
    recovery_probability = 0.5
    outbreak_probability = 0
    N_t = 1000

    #seed = 10010
    random.seed(seed)

    G = mhr_graph(B,L,k=7,xi=0.4)
    T = mhr_tree_graph(B,L,k=7,xi=0.4)

    k1 = p*(B-1)*(1+(xi-xi**L)/(1.-xi))
    k = mean(list(G.degree().values()))
    print(k1,k)
    i0 = 0.01
    N_i0 = int(i0 * N)
    patients_zero = random.choice(G.nodes(),N_i0,replace=False)

    simulation =  SIR(G,infection_probability,recovery_probability,outbreak_probability,save_node_changes=True,patients_zero=patients_zero)
    simulation2 = SIS(G,infection_probability,recovery_probability,outbreak_probability,save_node_changes=True,patients_zero=patients_zero)
    simulation.step(N_t=N_t)

    N_t = simulation.T+1
    t = arange(N_t)

    simulation2.step(N_t=simulation.T*3)

    fig2, (ax3,ax4) = pl.subplots(1,2,figsize=(20,10))
    fig, (ax1,ax2) = pl.subplots(1,2,figsize=(20,10))
    points, lines = draw_hierarchical_edge_bundling(G,T,B,L,ax2,draw_tree=False,colors=ones((int(N),3)))

    S = array(simulation.S)/N
    I = array(simulation.I)/N
    R = 1-array(simulation.S)/N-array(simulation.I)/N
    Sp, = ax1.plot([],[],'k')
    Ip, = ax1.plot([],[],'r')
    Rp, = ax1.plot([],[],'g')
    i_star, = ax1.plot([0,N_t],2*[1-recovery_probability/infection_probability/k],'c')
    ax1.set_xlim([0,N_t])
    ax1.set_ylim([0,1])
    
    pos = graphviz_layout(G,prog='neato')
    #pos = nx.circular_layout(G,prog='neato')
    nx.draw_networkx(G,pos,with_labels=False,ax=ax3,node_size=2)
    deg = list(G.degree().values())
    ax4.hist(list(G.degree().values()),bins=arange(min(unique(deg)),max(unique(deg))+1))


    def init():
        for n in simulation.newly_infected[0]:
            points[n].set_mfc('r')
        for n in simulation.newly_recovered[0]:
            points[n].set_mfc('g')        
        for n in simulation.newly_susceptible[0]:
            points[n].set_mfc('w')        
        for e in simulation.activated_edges[0]:
            lines[e].set_c('r')            
        Sp.set_data([],[])
        Ip.set_data([],[])
        Rp.set_data([],[])
        return tuple(points.values())+tuple(lines.values())+(Sp,Ip,Rp,)

    def animate(i):
        for n in simulation.newly_infected[i]:
            points[n].set_mfc('r')
        for n in simulation.newly_recovered[i]:
            points[n].set_mfc('g')        
        for n in simulation.newly_susceptible[i]:
            points[n].set_mfc('w')        
        for e in simulation.activated_edges[i]:
            lines[e].set_c('r')            
        if i>3:
            for e in simulation.activated_edges[i-4]:
                lines[e].set_c('k')            

        Sp.set_data(t[:i],S[:i])        
        Ip.set_data(t[:i],I[:i])        
        Rp.set_data(t[:i],R[:i])        
        return tuple(points.values())+tuple(lines.values())+(Sp,Ip,Rp,)

            
    anim = animation.FuncAnimation(fig,animate,frames=len(t)-1,interval=1,init_func=init,blit=False)
    ndx = get_node_ids_for_module(3,0,L,B)
    Sp.set_data(t[:N_t],S[:N_t])        
    Ip.set_data(t[:N_t],I[:N_t])        
    Rp.set_data(t[:N_t],R[:N_t])        

    t_shifted = (t[:-1]+t[1:])/2.    
    dS = diff(S)
    dI = diff(I)
    dR = diff(R)

    fig, (ax5,ax6) = pl.subplots(1,2)
    ax5.plot(t_shifted,dS,'k')
    ax5.plot(t_shifted,dI,'r')
    ax5.plot(t_shifted,dR,'g')


    S = array(simulation2.S)/N
    I = array(simulation2.I)/N
    ax6.plot(arange(simulation2.T+1),S)
    ax6.plot(arange(simulation2.T+1),I)
    #anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])        
    pl.show()
