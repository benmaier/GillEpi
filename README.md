BEWARE: THIS PACKAGE IS DEPRECATED

Please use [tacoma](github.com/benmaier/tacoma) instead, which is maintained and works better in every way.

# GillEpi

Provides pure Python classes to simulate epidemics on (potentially time varying) networks using a Gillespie stochastic simulation algorithm or the standard ABM SIS, SIR models.

## Install 

### Development version

    $ make

### Standard

    $ make install

## Examples

### List of classes

```python
from GillEpi import SI, SIS, SIR, SIRS
from GillEpi.agent_based_epidemics import SIR as AB_SIR
from GillEpi.agent_based_epidemics import SIS as AB_SIS
```

Find out the functionality using Python's help function and the examples below.

### Standard

```python
import GillEpi
import matplotlib.pyplot as pl
import networkx as nx

N = 100
k = 8
p = k / (N-1.0)
G = nx.fast_gnp_random_graph(N, p)

R0 = 1.5
recovery_rate = 1.0
infection_rate = R0 * recovery_rate / k
tmax = 1000

sis = GillEpi.SIS(
                  G,
                  infection_rate = infection_rate,
                  recovery_rate = recovery_rate,
                 )

# simulate
sis.simulate(tmax)

# plot infected cluster
i, t = sis.get_i_of_t()
pl.step(t,i)

pl.show()
```

### Agent-based model

I'm not a big fan of the node-centric ABM since the reaction `S+I - > I+I` is not being reflected with the right rates.

```python
from GillEpi.agent_based_epidemics import SIS as ABM_SIS
import matplotlib.pyplot as pl
import numpy as np
import networkx as nx

N = 100
k = 8
p = k / (N-1.0)
G = nx.fast_gnp_random_graph(N, p)

R0 = 1.5
recovery_probability = 0.01
infection_probability = R0 * recovery_probability / k
tmax = 1000

sis = ABM_SIS(
              G,
              infection_probability = infection_probability,
              recovery_probability = recovery_probability,
              patients_zero = [0,1,2,45,34],
             )

# simulate
sis.step(tmax)

# plot infected cluster
i = sis.I
t = sis.time

pl.step(t, i)

pl.show()
```

### Dynamic network

As an example for time-varying networks, I use the flockwork model (https://github.com/benmaier/flockworks).

```python
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
```
