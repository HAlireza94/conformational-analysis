import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import pandas as pd
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from lmfit import Model
from sympy import symbols, solve
import openpyxl 
from sklearn.metrics import mean_absolute_error,mean_squared_error
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import optimize
import sympy as sp
import math as m
from matplotlib import pyplot as plt, patches

import numpy as np
import math
from ase.io import read, write
from ase.geometry import get_dihedrals
from ase import Atoms

molecule = read('butane.gro', format='gromacs')
forcefield = []
with open ('butane.itp') as f:
    for line in f:
        p = line.split()

        forcefield.append(p)
atomNumber = molecule.get_global_number_of_atoms()


atomtypes, moleculetype, bonds, angles, dihedrals, pairs = [], [], [], [], [], []


for i in range(len(forcefield)):
    if len(forcefield[i]) == 3 and forcefield[i][1] == 'atomtypes':
        for j in range(i+1,i+atomNumber+1):
            # print(i)
            atomtypes.append(forcefield[j])
        print('================== Read AtomTypes ==================')        
    if len((forcefield[i])) == 3 and forcefield[i][1] == 'moleculetype':
        for k in range(i+5,i+5+atomNumber):
            moleculetype.append(forcefield[k])
        print('================== Read MoleculeType ==================')

    if len((forcefield[i])) == 3 and forcefield[i][1] == 'bonds':
        for s in range(i+1,i+atomNumber):
            bonds.append(forcefield[s])
        print('================== Read Bonds ==================')

    if len((forcefield[i])) == 3 and forcefield[i][1] == 'angles':
        for s in range(i+2,i+2*atomNumber-2):
            angles.append(forcefield[s])
        print('================== Read Angles ==================')

    if len((forcefield[i])) == 4 and forcefield[i][1] == 'PROPER':
        for s in range(i+2,i+2*atomNumber+1):
            dihedrals.append(forcefield[s])
        print('================== Read Dihedrals ==================')

    if len((forcefield[i])) == 3 and forcefield[i][1] == 'pairs':
        for s in range(i+1,i+2*atomNumber):
            pairs.append(forcefield[s])
        print('================== Read Pairs ==================')

def RB(psy, A0, A1, A2, A3, A4, A5):
    ANG = psy - 180
    RB = A0 + A1 * np.cos(np.deg2rad(ANG)) + A2 * (np.cos(np.deg2rad(ANG)))**2 + A3 * (np.cos(np.deg2rad(ANG)))**3 + A4 * (np.cos(np.deg2rad(ANG)))**4 + A5 * (np.cos(np.deg2rad(ANG)))**5
    return RB

angle = np.linspace(0,180,31)
energy_Q,energy_LJ, energy_RB, energy_bonds,energy_angle = [], [], [], [],[]
for j in range(len(angle)):
    
    molecule = read('butane.gro', format='gromacs')
    # molecule.set_pbc([False, False, False])  # Enable PBC in all directions (x, y, z)
    # molecule.wrap()
    bx = by = bz = molecule.get_cell()[0][0] * 0.1
    molecule.set_dihedral(3,2,1,0, angle=angle[j], mask=[1,1,0,0,1,1,1,1,1,0,0,0,0,0]) #mask=[1,1,0,0,1,1,1,1,1,0,0,0,0,0]
    write('butane-'+str(j)+'.gro', molecule)
    rb = 0
    for i in range(len(dihedrals)):
        phi = molecule.get_dihedral(int(dihedrals[i][0])-1, int(dihedrals[i][1])-1, int(dihedrals[i][2])-1, int(dihedrals[i][3])-1, mic=False)
        c0, c1, c2, c3 = float(dihedrals[i][5]), float(dihedrals[i][6]), float(dihedrals[i][7]), float(dihedrals[i][8])
        c4, c5 = 0, 0
        rb = rb + RB(phi,c0,c1,c2,c3,c4,c5)
    # print(rb)
    energy_RB.append(rb)

    U_b = 0
    for ii in range(len(bonds)):
        U_b = U_b + 0.5*float(bonds[ii][4])*(float(bonds[ii][3])-0.1*molecule.get_distance(int(bonds[ii][0])-1,int(bonds[ii][1])-1))**2
    
    energy_bonds.append(U_b)


    U_ang = 0
    for ik in range(len(angles)):
        U_ang = U_ang + 0.5*float(angles[ik][5])*(np.deg2rad(float(angles[ik][4])-molecule.get_angle(int(angles[ik][0])-1,int(angles[ik][1])-1, int(angles[ik][2])-1)))**2
    energy_angle.append(U_ang)



    epsilon, sigma = [], []
    for sd in range(len(moleculetype)):
        for sf in range(len(atomtypes)):
            if atomtypes[sf][0] == moleculetype[sd][1]:
            
                epsilon.append(float(atomtypes[sf][6]))
                sigma.append(float(atomtypes[sf][5]))

    
    # molecule.set_pbc([True, True, True])  # Enable PBC in all directions (x, y, z)
    # molecule.wrap()
    bx = by = bz = molecule.get_cell()[0][0] * 0.1

    x = molecule.get_positions()[:,0] * 0.1 # all nm
    y = molecule.get_positions()[:,1] * 0.1
    z = molecule.get_positions()[:,2] * 0.1 

    rc2 = 1.3**2 # nm^2
    numberPairs_14 = atomNumber*(atomNumber+1)/2 # all possible combinations


    en = 0
    count = 0
    for i in range(len(pairs)):
    
    
        dx = x[int(pairs[i][0])-1] - x[int(pairs[i][1])-1]
        dy = y[int(pairs[i][0])-1] - y[int(pairs[i][1])-1]
        dz = z[int(pairs[i][0])-1] - z[int(pairs[i][1])-1]
    
    
        DX = dx - bx * np.round(dx/bx)
        DY = dy - by * np.round(dy/by)
        DZ = dz - bz * np.round(dz/bz)

            
        r2 = DX**2 + DY**2 + DZ**2
    
    
        if r2 < rc2:
            count += 1
            epsilon_ij = np.sqrt(epsilon[int(pairs[i][0])-1] * epsilon[int(pairs[i][1])-1])
            sigma_ij = 0.5 * ( sigma[int(pairs[i][0])-1] + sigma[int(pairs[i][1])-1])                
            C6 = (sigma_ij**2/r2)**3            
            C12 = C6**2
            lj_rc = 4 * epsilon_ij * ((sigma_ij**2/rc2)**6 - (sigma_ij**2/rc2)**3)
            LJ = 4 * epsilon_ij * (C12 - C6)
            LJ -= lj_rc * C6
                
            en += LJ
                
            
    energy_LJ.append(en/2)
    charges = []
    for i in range(len(moleculetype)):
        charges.append(float(moleculetype[i][6]))
    
    en = 0
    for i in range(len(pairs)):
    
    
        dx = x[int(pairs[i][0])-1] - x[int(pairs[i][1])-1]
        dy = y[int(pairs[i][0])-1] - y[int(pairs[i][1])-1]
        dz = z[int(pairs[i][0])-1] - z[int(pairs[i][1])-1]
         
    
        DX = dx - bx * np.round(dx/bx)
        DY = dy - by * np.round(dy/by)
        DZ = dz - bz * np.round(dz/bz)

            
        r = np.sqrt(DX**2 + DY**2 + DZ**2)
        r2 = r**2

        qij = charges[int(pairs[i][0])-1] * charges[int(pairs[i][1])-1]

        f = 138.935458 # 1/(4*pi()*epsilon0)
        er = 1

        if r2 < rc2:
            coulomb = f*qij/(er*r)   
            en += coulomb
    energy_Q.append(en/2)

potential = []
for i in range(len(energy_angle)):
    potential.append(energy_angle[i] + energy_bonds[i] + energy_LJ[i] + energy_RB[i] + energy_Q[i])

pot = []
for i in range(len(potential)):
    pot.append((potential[i] - np.min(potential))/4.184)
angle1 = angle + 180
pot1 = pot[::-1]

plt.figure(dpi=250)
plt.plot(angle,pot,'-.k')
plt.plot(angle1,pot1,'-.k')
plt.xlim([0,360])
plt.ylim([0,10])
plt.title('OPLS-AA from LigPargen')
plt.xlabel('$C0-C01-C02-C03$ $(\degree)$')
plt.ylabel('$Relative$ $Energy$ $(kcal/mol)$')
plt.savefig('Potential Engergy as a function of Phi')


