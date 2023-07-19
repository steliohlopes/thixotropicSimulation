from dolfin import *


# Deformation Tensor
def DD(u):
    # Cartesian
    D = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    return D


# Stress Tensor
def TT(u, p, mu):
    # Cartesian
    T = 2 * mu * DD(u) - p * Identity(len(u))
    return T


def gammaDot(u):
    return pow(2 * inner(DD(u), DD(u)), 0.5)


# def eta(k,nPow,u):
#     eps=DOLFIN_EPS
#     return k*pow(gammaDot(u)+eps,nPow-1)
