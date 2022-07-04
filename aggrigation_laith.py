# creating  a function to calculate no. of ode required
import pandas as pd
from typing import Union
from sympy import summation, IndexedBase, Sum
import sympy as sympy

from sympy import symbols, Function, Eq, diff, Basic, Atom, IndexedBase, var, dsolve
import math
from pandas import DataFrame

dp = 1 * 10 ** (-9)
da = 1 * 10 ** (-5)
df: int = 3


def rhs(da, dp, df):
    return (da / dp) ** df;


na = int(rhs(da, dp, df))


# print("na=", na)


def nodes(n):
    return pow(2, math.ceil(math.log(n) / math.log(2)))


m = 0
for j in range(0, 200):
    if 2 ** j == nodes(na):
        m = j
        print("no.of odes to be formulated are:", m)
# radius list
# lets create some function
R = 2
d1 = dp
eps = 0.6366


def di(d1_arg: object, R_arg: int, i_arg: int, eps_arg: float) -> object:
    dii = (((R_arg ** (i_arg - 1)) / eps_arg) ** (1 / 3)) * d1_arg  # diameter of ith aggregate
    return dii


rad_d = []
# check-point
N = m + 3

for i in range(0, N):
    # print(i)
    radd = di(d1, R, i, eps)
    rad_d.append(radd)
    i = i + 1
# 0 - 42 actual
# 2 - 40 expected
# print(rad_d)
# len(rad_d)

# stage 2
di = rad_d
dj = rad_d
# check-point
d_combined = []
# d_combined = [[a, b] for a in di for b in dj]
for a in di:
    for b in dj:
        d_combined.append([a, b])

# print(len(d_combined))


# Collision Kernel, constants values for K function
Rg = 8.3 * (10 ** -3)
T = 298
B: float = 0.0014
um = 0.16


# Collision Kernel Function
def k(dii, djj):
    ##generating symbols for concentration variables of all aggregates
    # numSpecies = m + 1
    kij: int = (2 * Rg * T * B / (3 * um)) * (((dii + djj) ** 2) / (dii * djj))
    return kij


# if made a loop through d_combined then the result will be same
K_kernel = []
for i in rad_d:
    for j in rad_d:
        K_kernel.append(k(i, j))

K_kernel = pd.DataFrame(K_kernel, columns=["k"])
df = DataFrame(d_combined, columns=["di", "dj"])
df["k"] = K_kernel
# df
# print(df)

dff = df["k"]
# print(dff)
# ddf = pd.DataFrame(ddf[0;40]).T
K = []
# print(len(K_kernel))
# [ [0-42],[43,84] ... ]
for i in range(0, N):
    # startIndex = i * 42
    # qf = list(dff[startIndex:i * 40 + 40])
    # print(i)
    spacer = 1
    if i == 0:
        spacer = 0
    startIndex = i * 42 + spacer  # 0 , 43 , 85
    endIndex = 42 * i + 42 + spacer  # 42, 85 , 127
    # print([startIndex, endIndex])
    qf = list(dff[startIndex:endIndex])
    K.append(qf)

kk = K #pd.DataFrame(K)
print(kk)

# stage3
# generating symbols for concentration variables of all aggregates

numSpecies = N
# print(numSpecies)
C = [sympy.symbols("N%d" % numSpecies, cls=Function, Function=True) for numSpecies in range(40)]
# print(C)
t = symbols("t")

# N = sympy.symbols("N%d" % numSpecies)

eq = []


def K_Gm1(i, j):
    return kk[i][j]


for i in range(1, len(C)):

    if i == 1:
        eq1 = Eq(diff(C[i](t), t), (K_Gm1(i - 1, i - 1) / 2) * (10 ** 2))
        eq.append(eq1)
    else:
        eq1 = Eq(diff(C[i](t), t), (K_Gm1(i - 1, i - 1) / 2) * (C[i - 1](t) ** 2))
        eq.append(eq1)


print(eq[1:5])

# soln = dsolve((eq[ 1:5 ]), hint="default")
# soln
# print(soln)
# ------------------------------------------------------------------------------------------

# stage4

numSpecies = N
C = []
for num_Species in range(0, N):
    C.append(sympy.symbols("N%d" % num_Species, integer=True))
t = symbols("t")


def K_Gm2(i, j):
    return kk[i][j]

# i 1:5
# j 1:5
#1,1 1,2 1,3 1,4 1,5
#2,1 2,2 2,3 2,4 2,5
#3,1 3,2 3,3 3,4 3,5
#4,1 4,2 4,3 4,4 4,5
#5,1 5,2 5,3 5,4 5,5
#j = 1  # this code could change everything
eq2 = []
J = sympy.Symbol('j')
I= sympy.symbols('i')

def symSum(lowerLimit:int,upperLimit:int,iValue:int):
    counter = 0
    for j in range(lowerLimit,upperLimit+1):
        counter = counter + K_Gm1(iValue - 1, j) * C[j] * t * ((2 ** (j - 1)) / (2 ** (iValue - 1) - (2 ** (iValue - 2))))
    return counter

for i in range(3, N):
    eq_index1 = Eq(
        diff(Function(C[i])(t), t), # i:x
        C[i - 1] *
        t *
        symSum(1,i-2,i)
        # summation(K_Gm1(i - 1, j) * (C[j] * t * ((2 ** (J - 1)) / (2 ** (i - 1) - (2 ** (i - 2))))), ( J ,1, i))
    )
    eq2.append(eq_index1)
    if i == 5:
        print(eq_index1)

# (C[j], C[2], C[i-3]) mean C OF FIRST LIMIT AND SECOND NOT BETWEEN
