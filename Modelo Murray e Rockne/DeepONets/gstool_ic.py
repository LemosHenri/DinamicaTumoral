import gstools as gs
import numpy as np 
import matplotlib.pyplot as plt

from numpy.random import randint
from gstools.random import MasterRNG

L = 20
N_X = 100
x = np.linspace(0, 1, N_X)


branch = L**3 * np.exp(-100 * x**2) 
points = randint(100, size = 10)

samples = [[x[0], branch[0]]]
for i in points: samples.append([x[i], branch[i]])
samples = np.array(samples)

seed = MasterRNG(123123)
def Random_functions(x, samples):

    points = randint(100, size = 5)

    samples = [[x[0], branch[0]]]
    for i in points: samples.append([x[i], branch[i]])
    samples = np.array(samples)
    model = gs.Gaussian(dim=1, var=10, len_scale=10)
    krige = gs.krige.Simple(model, samples[:, 0], samples[:, 1])
    srf = gs.CondSRF(krige, seed = seed())
    return srf.structured([x])

N_ICS = 250
ICS = np.zeros((N_ICS, N_X))
for fn in range(N_ICS):
    x_values = Random_functions(x, samples)
    ICS[fn, :] = x_values

plt.plot(x, ICS[10])
plt.plot(x, ICS[20])
plt.plot(x, ICS[30])
plt.plot(x, ICS[40])
plt.plot(x, ICS[50])
#plt.scatter(samples[:, 0], samples[:, 1])
plt.grid()
plt.show()