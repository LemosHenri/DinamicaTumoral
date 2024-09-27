import numpy as np

def Params(D: float, p: float, L: float, tf: int, c0: float, a: float, ab: float) -> dict:

    Dl = D / (p * L**2)

    return {'D': D, 'p': p, 'L': L, 'c0': c0, 'Dl': Dl, 'alpha': a, 'ab': ab, 'tf': tf}

def R(d: float, apply: bool, params: dict) -> float:
    if apply:
        return -1. + (1. - np.exp(-params['alpha'] * (d + (d**2 / params['ab'])))) / params['p']
    else:
        return -1.
    
def Matrix(dose: float, apply: bool, params: dict, dx: float, dt: float, size: int):

    LAMBDA = (params['D'] * dt) / dx**2
    GAMMA_1 = 1 + (R(dose, apply, params) * dt) / 2 + LAMBDA 
    GAMMA_2 = 1 - (R(dose, apply, params) * dt) / 2 - LAMBDA 

    A = np.zeros((size, size))
    B = np.zeros((size, size))

    for i in range(1, size-1):
        A[i,i] = GAMMA_1
        A[i,i+1] = -LAMBDA / 2
        A[i,i-1] = -LAMBDA / 2

        B[i,i]  = GAMMA_2
        B[i,i+1] = LAMBDA / 2
        B[i,i-1] = LAMBDA / 2


    A[0, 0] = - 3 / (4 * dx)
    A[0, 1] = 1 / dx
    A[0, 2] = -1 / (4 * dx)

    A[-1, -1] = 3 / (4 * dx)
    A[-1, -2] = -1 / dx
    A[-1, -3] = 1 / (4 * dx)

    B[0, 0] = 3 / (4 * dx)
    B[0, 1] = -1 / dx
    B[0, 2] = 1 / (4 * dx)

    B[-1, -1] = -3 / (4 * dx)
    B[-1, -2] = 1 / dx
    B[-1, -3] = -1 / (4 * dx)

    return A, B

def Solve(days, doses, params, size_x, size_t, step = 1):

    X = np.linspace(0, 1, size_x) 
    T = np.linspace(0, params['tf'] * params['p'], size_t) 
    X_K = params['L']**3 * np.exp(-100 * X**2)
    S = np.zeros((size_t, size_x))

    for it in range(size_t):

        if it in days:
            index = days.index(it)
            A, B = Matrix(doses[index], True, params, X[1], T[1], size_x)

        else: 
            A, B = Matrix(0., False, params, X[1], T[1], size_x)

        S[it] = X_K
        E = B @ X_K
        X_K = np.linalg.solve(A, E)
        if it % step == 0:
            print(f'Iteração: {it}/{size_t}')

    return X, T, S

def Tumor_radius(s, x, size_x, size_t, params):
    
    R = np.zeros(size_t)

    for j in range(size_t):
        
        for i in range(size_x):

            if s[j, i] <= params['c0'] * 0.6126 and i > 3:
                R[j] = x[i] * params['L']
                break
    return R
