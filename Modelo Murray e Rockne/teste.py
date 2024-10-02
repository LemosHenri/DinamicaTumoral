import numpy as np

D = 1
dt = 2
dx = 3

LAMBDA = -4 * (D * dt) / dx**2 + dt + 1
GAMMA = D * dt / dx**2

N_S = 3
M = np.zeros((N_S**2, N_S**2))

M[0,0] = 6.
M[0,1] = -4.
M[0,2] = 1.
M[0,N_S] = -4.
M[0,2*N_S] = 1.

M[-1,-1] = 6.
M[-1,-2] = -4.
M[-1,-3] = 1.
M[-1,-N_S-1] = -4.
M[-1,-2*N_S-1] = 1.


for i in range(1, N_S**2 - 1):

        if i< N_S:
                M[i, i] = - 3 / (4 * dx)
                M[i, i + N_S] = 1 / dx
                M[i, i + 2*N_S] = -1 / (4 * dx)
                
        if i > N_S:
                M[i,i] = LAMBDA
                M[i,i+1] = GAMMA
                M[i,i-1] = GAMMA
                #M[i,i+N_S] = GAMMA
                #M[i,i-N_S] = GAMMA

        if i > N_S**2 - N_S:
                M[N_S - i, i] = 3 / (4 * dx)
                M[N_S - i, i - N_S] = - 1 / dx
                M[N_S - i, i - 2*N_S] = 1 / (4 * dx)

print(M)