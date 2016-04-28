from flockworks import flockwork
import GillEpi
import pylab as pl

# initialize time-varying network
F = flockwork(Q=0.7,N=100,k0=2)
F.equilibrate()

#initialize SIR simulation
sir = GillEpi.SIR(
                  F.G,
                  infection_rate = 1.,
                  recovery_rate = 1.,
                  rewiring_rate = 1.,
                  rewire_function = F.rewire,
                  mean_degree_function = F.mean_degree
                 )

# simulate
sir.simulate()

# initialize analysis
fig, ax = pl.subplots(2,1)

# plot susceptible cluster
s, t = sir.get_s_of_t()
ax[0].step(t,s)

# plot resistant cluster
r, t = sir.get_r_of_t()
ax[0].step(t,r)

# plot infected cluster
i, t = sir.get_i_of_t()
ax[0].step(t,i)

# plot basic reproduction number
R0,t = sir.get_R0_of_t()
ax[1].step(t,R0)

pl.show()
