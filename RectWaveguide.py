import cmath

import scipy.constants
import sympy as sy
from matplotlib import pyplot as plt, figure
from itertools import repeat
import numpy as np
from sympy import symbols, solve, I, nsolve, solveset, S, Eq, pprint
from scipy import constants as cc
import cmath
# import psutil
#
#
# def checkIfProcessRunning(processName):
#     # Iterate over the all the running process
#     for proc in psutil.process_iter():
#         try:
#             # Check if process name contains the given name string.
#             if processName.lower() in proc.name().lower():
#                 return True
#         except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
#             pass
#     return False
#
#
# def findProcessIdByName(processName):
#     listOfProcessObjects = []
#     # Iterate over the all the running process
#     for proc in psutil.process_iter():
#         try:
#             pinfo = proc.as_dict(attrs=['pid', 'name', 'create_time'])
#             # Check if process name contains the given name string.
#             if processName.lower() in pinfo['name'].lower():
#                 listOfProcessObjects.append(pinfo)
#         except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
#             pass
#     return listOfProcessObjects
#
#
# def is_port_in_use(port):
#     import socket
#     with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
#         return s.connect_ex(('localhost', port)) == 0
#
#
# def scan_ports(x1, x2):
#     for i in range(x1, x2):
#         k = is_port_in_use(i)
#         print(k)
#         print(i)
#         if k:
#             return i
#         else:
#             pass
#
#
# def showProcAndPorts(i):
#     import os
#     if i == 0:
#         a = os.popen('wmic process get description, processid').read()
#     else:
#         a = os.popen('netstat -nb').read()
#     print("\n Connections ", a)
#
#
# Port = 0
# listProc = []
# if checkIfProcessRunning('fimmwave'):
#     print('Yes a fimmwave process was running')
#     listOfProcessIds = findProcessIdByName('fimmwave')
#     if len(listOfProcessIds) > 0:
#         print('Process Exists | PID and other details are')
#         for elem in listOfProcessIds:
#             processID = elem['pid']
#             listProc.append(processID)
#             processName = elem['name']
#             print((processID, processName))
#         Port = min(listProc)
#     else:
#         print('No Running Process found with given text')
#
# else:
#     print('No fimmwave process was running')

#showProcAndPorts(0)

## Reactangle Waveguides - 4.25um:
def power(list, pow):
    return [x ** pow for x in list]


def ZeroDiv(x, y):
    try:
        return x / y
    except ZeroDivisionError:
        return 0


def modeColor(var):
    if var == 0:
        color = 'blue'
        return color
    if var == 1:
        color = 'red'
        return color
    if var == 2:
        color = 'green'
        return color
    if var == 3:
        color = 'purple'
        return color
    if var == 4:
        color = 'orange'
        return color


def dNeff(d, m, p):
    first_part = -(k0 * d * sy.sqrt(n1 ** 2 - neffm ** 2))
    second_part = m * sy.pi
    b = sy.atan(((n1 / nb) ** (2 * p)) * sy.sqrt((neffm ** 2 - nb ** 2) / (n1 ** 2 - neffm ** 2)))
    c = sy.atan(((n1 / nc) ** (2 * p)) * sy.sqrt((neffm ** 2 - nc ** 2) / (n1 ** 2 - neffm ** 2)))
    last_part = sum([b, c])
    expr = first_part + second_part + last_part
    return expr


def wNeff(w, neff1, neff2, r, p):
    first_part = k0 * w * sy.sqrt(n1 ** 2 - Nmr ** 2)
    second_part = r * sy.pi
    last_part = 2 * sy.atan((neff2[0] / neff1[0]) ** (2 * p)) \
                * sy.sqrt((Nmr ** 2 - neff1[0] ** 2) / (neff2[0] ** 2 - Nmr ** 2))
    expr = -first_part + second_part + last_part
    return expr


nAir = nc = n0 = n3 = 1.0003
nGe = ng = n1 = nf = 3.979
nSi = ns = nb = n2 = 3.425
wavelength = 4250  # [nm]
w = 1000.  # [nm]
hng = 1000.  # [nm]
hns = 2000.  # [nm]
t = 500  # [nm]
d = hng - t  # [nm]
neffm = symbols('neffm', real=True)
Nmr = symbols('Nmr', real=True)
# Effective index method: ASYMMETRIC SLAB WAVEGUIDE DESIGN
# Optical rib waveguides based on sol-gel derived silica–titania films
k0 = 2 * sy.pi / wavelength
pTE = 0
pTM = 1
p = pTE
i = 0
j = 0
cut_off = []
k = 4500
o = 10
dmin_modeTE = []
dmin_modeTM = []
index_modeTE = []
index_modeTM = []
sol_mode = []
neff1 = []
neff2 = []
fig1, ax1 = plt.subplots()
for p in [pTE, pTM]:
    for m in range(2, -1, -1):
        verr = 0
        color = modeColor(m)
        x = []
        y = []
        sol = []
        for d in range(k, 0, -o):
            expr = dNeff(d, m, p)
            x.append(d)
            try:
                sol.append(nsolve(expr, neffm, 3.7, real=True).as_real_imag()[0])
            except ValueError:
                if verr == 0 and p == pTE:
                    dmin_modeTE.append(x[i])
                    index_modeTE.append(j)
                    verr += 1
                else:
                    pass
                if verr == 0 and p == pTM:
                    dmin_modeTM.append(x[i])
                    index_modeTM.append(j)
                    verr += 1
                else:
                    pass
                sol.append(sol[i - 1])
            sol_mode.append(sol[i])
            y.append(sol[i])
            if len(dmin_modeTE) <= 2:
                pass
            else:
                if d == dmin_modeTE[1] and m == 1:
                    neff1.append(sol_mode[j])
            i += 1
            j += 1
        if p == 0:
            mode = 'TE' + str(m)
            plt.plot(x, y, '-', color=color, label=mode)
        if p == 1:
            mode = 'TM' + str(m)
            plt.plot(x, y, '--', color=color, label=mode)
        i = 0

neff2.append(sol_mode[index_modeTE[2]])
print('neff2',neff2)
cut_off.extend(repeat(min(y), int(k / o)))
plt.plot(x, cut_off, '--', color='gold', label='cutoff')
plt.xlabel('h [nm]')
plt.ylabel('Neff')
plt.grid()
plt.legend()
print("Cut_Off height for TE mode{}".format(dmin_modeTE))
print("Cut_Off height for TM mode{}".format(dmin_modeTM))
print("Cut_Off index for TE mode{}".format(index_modeTE))
print("Cut_Off index for TM mode{}".format(index_modeTM))
for xc in dmin_modeTE:
    plt.axvline(x=xc, linestyle='--', color='red', alpha=0.3)
for xc in dmin_modeTM:
    plt.axvline(x=xc, linestyle='--', color='blue', alpha=0.3)
i = 0
j = 0
sol_1 = []
sol_2 = []

# RIB WAVEGUIDE DESIGN : WIDTH
plt.figure(2)
for p in [pTE, pTM]:
    for r in range(0, 3):
        verr = 0
        color = modeColor(r)
        x = []
        y = []
        sol = []
        for w in range(4000, 100, -100):
            expr = wNeff(w, neff2, neff1, r, p)
            x.append(w)
            try:
                sol.append(nsolve(expr, Nmr, 3.7, real=True).as_real_imag()[0])
            except ValueError:
                sol.append(abs(sol[i - 1]))
            y.append(abs(sol[i]))
            i += 1
        if p == 0:
            plt.plot(x, y, '-', color=color)
        if p == 1:
            plt.plot(x, y, '--', color=color)
        i = 0
plt.grid()
plt.xlabel('w [nm]')
plt.ylabel('Neff')
# WAVE PROPAGATION SLAB WG 2D
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
E0 = H0 = 1
hng = 2000
edge1 = k0 * ns  # max([nc, ns])  # ns is bigger
edge2 = k0 * nf
beta_arr = np.arange(edge1, edge2, 0.01 * edge1)
beta = np.mean(beta_arr)
print(beta_arr)
print(beta)
print(k0 * nf)
gamma_c = Y3 = cmath.sqrt(beta ** 2 - k0 ** 2 * nc ** 2)
gamma_s = Y2 = cmath.sqrt(beta ** 2 - k0 ** 2 * ns ** 2)
kfx = h1 = cmath.sqrt(k0 ** 2 * nf ** 2 - beta ** 2)
a = hng / 2
W = Y2 * a
U = h1 * a
Wprim = Y3 * a
Em0 = []
Em1 = []
Em2 = []
Em3 = []
Hm0 = []
Hm1 = []
Hm2 = []
Hm3 = []
X = []
u = (0 * cmath.pi / 2) + (1 / 2) * cmath.atan(W / h1 * a) + (1 / 2) * cmath.atan(Wprim / h1 * a)
for m in range(0, 4):
    for x in range(-30000, 30000, 1):
        x = x * 0.1
        if m == 0:
            X.append(x)
        u = (m * cmath.pi / 2) + (1 / 2) * cmath.atan(W / u) + (1 / 2) * cmath.atan(Wprim / u)
        # phiTE = cmath.atan(gamma_s / kfx)
        phiTE = (m * cmath.pi / 2) + (1 / 2) * cmath.atan(W / u) - (1 / 2) * cmath.atan(Wprim / u)
        # phiTM = cmath.atan(nf**2 * Y2 / (ns**2 * kfx))
        phiTM = (m * cmath.pi / 2) - (1 / 2) * cmath.atan(W / u) + (1 / 2) * cmath.atan(Wprim / u)
        if x > a:
            Ey1 = E0 * cmath.cos(kfx * a - phiTE) * cmath.exp(-Y3 * (x - a))  # *cmath.exp(-I*beta) <- complex form
            Hy1 = E0 * cmath.cos(kfx * a - phiTM) * cmath.exp(-Y3 * (x - a))  # *cmath.exp(-I*beta) <- complex form
            if m == 0:
                Em0.append(Ey1.real)
                Hm0.append(Hy1.real)
            if m == 1:
                Em1.append(Ey1.real)
                Hm1.append(Hy1.real)
            if m == 2:
                Em2.append(Ey1.real)
                Hm2.append(Hy1.real)
            if m == 3:
                Em3.append(Ey1.real)
                Hm3.append(Hy1.real)
        if -a <= x <= a:
            Ey2 = E0 * cmath.cos(kfx * x - phiTE)  # *cmath.exp(-I*beta) <- complex form
            Hy2 = E0 * cmath.cos(kfx * x - phiTM)
            if m == 0:
                Em0.append(Ey2.real)
                Hm0.append(Hy2.real)
            if m == 1:
                Em1.append(Ey2.real)
                Hm1.append(Hy2.real)
            if m == 2:
                Em2.append(Ey2.real)
                Hm2.append(Hy2.real)
            if m == 3:
                Em3.append(Ey2.real)
                Hm3.append(Hy2.real)
        if x < -a:
            Ey3 = E0 * cmath.cos(kfx * a + phiTE) * cmath.exp(Y2 * (x + a))  # *cmath.exp(-I*beta) <- complex form
            Hy3 = E0 * cmath.cos(kfx * a + phiTM) * cmath.exp(Y2 * (x + a))
            if m == 0:
                Em0.append(Ey3.real)
                Hm0.append(Hy3.real)
            if m == 1:
                Em1.append(Ey3.real)
                Hm1.append(Hy3.real)
            if m == 2:
                Em2.append(Ey3.real)
                Hm2.append(Hy3.real)
            if m == 3:
                Em3.append(Ey3.real)
                Hm3.append(Hy3.real)

ax1.plot(X, Em0, 'b-')
ax2.plot(X, Em1, 'r-')
ax3.plot(X, Em2, 'g-')
ax4.plot(X, Em3, '-', color='purple')
ax1.plot(X, Hm0, 'b--')
ax2.plot(X, Hm1, 'r--')
ax3.plot(X, Hm2, 'g--')
ax4.plot(X, Hm3, '--', color='purple')
ax1.axvline(x=a, linestyle='-.', color='black')
ax1.axvline(x=-a, linestyle='-.', color='black')
ax2.axvline(x=a, linestyle='-.', color='black')
ax2.axvline(x=-a, linestyle='-.', color='black')
ax3.axvline(x=a, linestyle='-.', color='black')
ax3.axvline(x=-a, linestyle='-.', color='black')
ax4.axvline(x=a, linestyle='-.', color='black')
ax4.axvline(x=-a, linestyle='-.', color='black')
plt.tight_layout()
ax1.grid(), ax2.grid(), ax3.grid(), ax4.grid()
plt.show()


## Effective index method: ASYMMETRIC RIB WAVEGUIDE DESIGN
# Large Single-Mode Rib Waveguides in GeSi-Si and Si-on-SO2
# +
# Dimension Dependent Silicon Photonics Rib Waveguide Characterisation
height1 = []
height2 = []
width1 = []
aE = []
bE = []
rE = []
aH = []
bH = []
rH = []
single_mode_a = []
single_mode_b = []
# plt.figure(3)
for j in range(0, 2):
    gammaHE = 1
    gammaEH0 = (n0 / n1) ** 2
    gammaEH2 = (n2 / n1) ** 2
    if j == 0:
        gamma0 = gammaHE
        gamma2 = gammaHE
    else:
        gamma0 = gammaEH0
        gamma2 = gammaEH2
    for i in range(0, 200, 1):
        a = i * 0.01
        b = 2 * 1 / 8.5
        r = 0.55  # 0.5 < r < 1
        h1 = 2 * b * r * wavelength
        h2 = 2 * b * wavelength
        width = 2 * a * wavelength
        # print(width)
        # print(h1)
        # print(h2)
        q = (gamma0 / sy.sqrt(n1 ** 2 - n0 ** 2)) + (gamma2 / sy.sqrt(n1 ** 2 - n2 ** 2))
        w1 = 4 * sy.pi * b / (q + 4 * sy.pi * b)
        w2 = 4 * sy.pi * r * b / (q + 4 * sy.pi * r * b)
        hi = 2 * b * wavelength / w1
        h0 = 2 * b * r * wavelength / w2
        delta1 = (hi / h0) ** 2 - 1
        delta2 = (w2 / r * w1) ** 2 - 1
        V = (sy.pi / 2) * (a * w1 / b) * sy.sqrt(delta2)
        Vs = (sy.pi / 2) * (1 + 0.3 * sy.sqrt(delta2))
        if j == 0:
            height1.append(h1)
            height2.append(h2)
            width1.append(width)
            aE.append(a)
            bE.append(b)
            rE.append(r)
        else:
            height1.append(h1)
            height2.append(h2)
            width1.append(width)
            aH.append(a)
            bH.append(b)
            rH.append(r)
        if ZeroDiv(a, b) <= 0.3 + r / (sy.sqrt(1 - r ** 2)):  # SINGLE MODE CONDITION (SMC) : Soref et al
            single_mode_a.append(a)
            single_mode_b.append(b)
        else:
            continue
        # Chan et al proposition for small rib wg -> W/H <= 0.05 + ((0.94 + 0.25*H)*r/sqrt(1-r**2)) - 2
        # 0 <= r <= 0.5 and 1 <= H <= 1.5um

aE_arr = np.array(aE, dtype=float)
bE_arr = np.array(bE, dtype=float)
aH_arr = np.array(aH, dtype=float)
bH_arr = np.array(bH, dtype=float)
list_divE = np.true_divide(aE_arr, bE_arr, out=np.zeros_like(aE_arr), where=bE_arr != 0, casting='unsafe')
list_divH = np.true_divide(aH_arr, bH_arr, out=np.zeros_like(aH_arr), where=bH_arr != 0, casting='unsafe')
# plt.plot(aE_arr, list_divE, 'r')
# plt.plot(aH_arr, list_divH, 'b--')
print(single_mode_a)
print(single_mode_b)
bidx = next((i for i, x in enumerate(single_mode_b) if x), None)
aidx = next((i for i, x in enumerate(single_mode_a) if x), None)
print(bidx)
# idx = np.nonzero(single_mode_b)
# bidx = idx
# print(bidx)
print('h1:slab height:', single_mode_b[bidx] * r * wavelength)
print('h2:waveguide height:', single_mode_b[bidx] * 2 * wavelength)
print('w1:waveguide width:', single_mode_a[-1] * wavelength)
# plt.grid()

# RIDGE WAVEGUIDE DESIGN
# The ridge waveguide (Figure 3.5) is a special case of a rib waveguide, although it may
# be also understood as a generalised rectangular waveguide. The ridge wave-guide  has  similar  single  mode
# properties  to  a  rectangular  waveguide. From  the effective index analysis, it follows that slightly larger
# single mode waveguide cross sections, than in the case of rectangular waveguides, can be achieved [3]. Ebeling,
# K.J., Integrated Optoelectronics: Waveguide Optics, Photonics, Semiconductors. 1993, Berlin: Springer Verlag

w_s = []
h_s = []
for w in range(0, 2000):
    scond_w = w * (4 * sy.pi / wavelength) * sy.sqrt(nf ** 2 - nc ** 2)
    if scond_w < 2:
        w_s.append(w)
    else:
        continue

for h in range(0, 2000):
    scond_h = h * (2 * sy.pi / wavelength) * sy.sqrt(nf ** 2 - nc ** 2)
    if scond_h < sy.pi / 2:
        h_s.append(h)
    else:
        continue
print(w_s)
print('W for square single mode', w_s[-1])
print(h_s)
print('H for square single mode', h_s[-1])

# Extension of Marcatili’s analytical approach forrectangular silicon optical waveguides
n4 = n5 = n3
m = 1
n = 1
d = w = 1000
b = a = hng = 2000
airh = 1000
freq = 70.539402 * 10E+12
# Xi = symbols('Xi', real=True)
# mi = symbols('mi', real=True)

step = 50
# Xi = 1
# mi = 1
x = 0
y = 0
plt.figure(4)
for Xi in range(-int(d), int(d), step):
    for mi in range(-int(b), int(b), step):
        omega = k0 / cc.c / nGe
        kx = m * sy.pi / d
        ky = n * sy.pi / b
        beta = sy.sqrt(n1 ** 2 * k0 ** 2 - kx ** 2 - ky ** 2)
        # print('Beta:', beta)
        gamma2 = sy.sqrt((n1 ** 2 - n2 ** 2) * k0 ** 2 - kx ** 2)
        gamma3 = sy.sqrt((n1 ** 2 - n3 ** 2) * k0 ** 2 - kx ** 2)
        gamma4 = sy.sqrt((n1 ** 2 - n4 ** 2) * k0 ** 2 - ky ** 2)
        gamma5 = sy.sqrt((n1 ** 2 - n4 ** 2) * k0 ** 2 - ky ** 2)

        # print(mI)
        A1 = 1000
        A2 = A1 * omega * cc.epsilon_0 * n1 ** 2 * ky / (beta * kx)

        # print('A2:', A2)
        A3 = A1 * sy.sin(kx * (Xi - d / 2))
        A4 = A2 * sy.cos(kx * (Xi - d / 2))
        A5 = A1 * sy.sin(kx * (Xi + d / 2))
        A6 = A2 * sy.cos(kx * (Xi + d / 2))
        A7 = A1 * sy.cos(ky * (mi - b / 2))
        # A7 = A1 * (1 + (k0 ** 2 * (n1 ** 2 - n4 ** 2)) / beta ** 2) * sy.cos(ky * (mi - b / 2))
        # print('A7:', A7)
        A8 = A2 * sy.sin(ky * (mi - b / 2))
        A9 = A1 * sy.cos(ky * (mi + b / 2))
        # A9 = A1 * (1 + (k0 ** 2 * (n1 ** 2 - n5 ** 2)) / beta ** 2) * sy.cos(ky * (mi + b / 2))
        A10 = A2 * sy.sin(kx * (Xi + b / 2))
        y = Xi
        if -d / 2 <= Xi <= d / 2:
            if -b / 2 <= mi <= b / 2:
                #Ez = A1 * sy.sin(kx * (x + Xi)) * sy.cos(ky * (y + mi))
                Ex = ((beta * A1 * kx + omega * cc.mu_0 * A2 * ky) / (n1 ** 2 * k0 ** 2 - beta ** 2)) * sy.cos(kx * (x + Xi)) * sy.cos(ky * (y + mi))
                plt.plot(Xi,abs(Ex), 'ro', alpha=0.2)
                #Ey = ((-beta * A1 * ky + omega * cc.mu_0 * A2 * kx) / (n1 ** 2 * k0 ** 2 - beta ** 2)) * sy.sin(kx * (x + Xi)) * sy.sin(ky * (y + mi))
                # plt.plot(mi, abs(Ey) * w - w / 2, 'bo', alpha=0.2)
                # plt.plot(mi,Ez,'bo', alpha=1)
        if -d / 2 <= Xi <= d / 2:
            if mi <= -b / 2:
                Ez = A7 * sy.sin(kx * (x + Xi)) * sy.exp(gamma4 * (y + b / 2))
                #plt.plot(mi, abs(Ez) - w / 2, 'go', alpha=1)
        if -d / 2 <= Xi <= d / 2:
            if mi >= b / 2:
                Ez = A9 * sy.exp(-gamma5 * (y - b / 2)) * sy.sin(kx * (x + Xi))
                #plt.plot(Xi, abs(Ez), 'bo', alpha=1)
        # plt.plot(y, A3, 'ro')
        # plt.plot(y, A2, 'ro')
        # plt.plot(y, A4, 'ro')
        # plt.plot(y, A5, 'ro')
        # plt.plot(y, A6, 'ro')
        # y = mi
        # plt.plot(y, A7, 'bo')
        # plt.plot(y, A8, 'bo')
        # plt.plot(y, A9, 'bo')
        # plt.plot(y, A10, 'bo')
        # A2 = A1 * beta * ky / omega * cc.mu_0 * kx
        # A3 = A1 * sy.sin(kx * (Xi - d / 2))
        # A4 = A2 * (1 + (k0 ** 2 * (n1 ** 2 - n2 ** 2)) / beta ** 2) * sy.cos(kx * (Xi - d / 2))
        # A5 = A1 * sy.sin(kx * (Xi + d / 2))
        # A6 = A2 * (1 + (k0 ** 2 * (n1 ** 2 - n3 ** 2)) / beta * 2) * sy.cos(kx * (Xi + d / 2))

plt.axvline(x=-d / 2, ymin=-hng, ymax=hng)
plt.axvline(x=d / 2, ymin=-hng, ymax=hng)
plt.axhline(y=-hng/2, xmin=-w, xmax=w)
plt.axhline(y=hng/2, xmin=-w, xmax=w)
plt.axhline(y=-hng/2 - hns, xmin=-w, xmax=w)
plt.axhline(y=hng/2 + airh, xmin=-w, xmax=w)
plt.grid()
plt.show()
