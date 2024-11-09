from PlateKirchhoff import fem_solver
from PlateKirchhoff.fe_system import fe_data, assemble_system
# from copy import deepcopy
import time
import numpy as np

import matplotlib.pyplot as plt
# from matplotlib import animation

import multidim_galerkin_pod.gen_pod_utils as gpu

import sys
# import os
sys.path.append('./PlateKirchhoff')  # source


# Finite element system
solverTic = time.time()

# coord = np.array([[0.,  0. ], [1.,  0. ], [1.,  0.7], [0.,  0.7]])
c_s = np.array([1., 0.7, 0.003])

# Structured nodal coordinates
nElem_x = 20
# Number of elements in x direction
nElem_y = 20
# Number of elements in y direction

x = np.linspace(0, c_s[0], nElem_x+1)
y = np.linspace(0, c_s[1], nElem_y+1)

coord = np.zeros(((nElem_y+1)*(nElem_x+1), 2))  # Two coords per node (x,y)
a_x = c_s[0] / nElem_x
a_y = c_s[1] / nElem_y
for m in range(nElem_x+1):
    for n in range(nElem_y+1):
        coord[m*(nElem_y+1)+n] = np.array([(m)*a_x, (n)*a_y])


# Generate element connectivity
# Four nodes per element (linear quadrilateral)
connect = np.zeros((nElem_y*nElem_x, 4), dtype=np.int32)
for m in range(nElem_x):
    for n in range(nElem_y):
        elem_connect = np.zeros(4, dtype=np.int32)
        elem_connect[0] = (m)*(nElem_y+1) + (n)
        elem_connect[1] = elem_connect[0] + (nElem_y+1)
        elem_connect[2] = elem_connect[1] + 1
        elem_connect[3] = elem_connect[0] + 1
        connect[m*nElem_y+n] = elem_connect

# Material properties
mat = np.array([7.0e+10, 3.4e-01, 2.7e+03])
material = np.zeros((nElem_x*nElem_y, 3))
# Emod, Poissons ratio, Density
# (Same material for all elements / homogeneous domain)
material[:, :] = mat

n_unkwn_elem_disp = connect.shape[0]*4
solverToc = time.time()
solverTime = solverToc - solverTic
print('>> Time for mesh gen: ' + str(round(solverTime, 3)) + ' seconds.')

# data structure for FEM
my_fe_data = fe_data()
my_fe_data.coord = coord
my_fe_data.material = material
my_fe_data.pDim = c_s
my_fe_data.connect = connect
my_fe_data.bc = np.array([])
my_fe_data.grid_shape = (nElem_x+1, nElem_y+1)

# Matrix assembly
solverTic = time.time()
# KMAT - stiffness matrix and MMAT - mass matrix
KMAT, MMAT = assemble_system(my_fe_data)
solverToc = time.time()
print('>> Time for assembling: ' + str(round(solverTime, 3)) + ' seconds.')
DMAT = np.zeros_like(MMAT)                # DMAT - damping matrix

# Load vector
n = KMAT.shape[0]
fload = np.zeros((n, 1))
fload[0, 0] = 1

# Frequency response
freq = np.linspace(10, 400, 391)

# defining the function


def output_func(_freq):
    yFEM = fem_solver.solve_dense_system_response(
        KMAT, MMAT, DMAT, fload, _freq)
    # filter z translation
    yFEM = np.real(yFEM[0::4, :])
    yFEM = yFEM.reshape(my_fe_data.grid_shape[0], my_fe_data.grid_shape[1], -1)
    return yFEM


# get snapshot tensor with physical coordinates (x,y,frequency)
# solverTic = time.time()
# snapshots_freq = output_func(freq)
# solverToc = time.time()
# print('>> Time for solving: ' + str(round(solverTime,3)) + ' seconds.')
# np.savez(f'./data/01_data.npz',snapshots_freq=snapshots_freq)

snapshot_data = np.load('./data/01_data.npz')
snapshots_freq = snapshot_data['snapshots_freq']

# plot the frequency response function at a plate node
plot_nodeid_x = 0
plot_nodeid_y = 0
plt.figure()
plt.semilogy(freq, np.abs(
    snapshots_freq[plot_nodeid_x, plot_nodeid_y, :]), linewidth=2)
plot_nodeid_x = -1
plot_nodeid_y = -1
plt.semilogy(freq, np.abs(
    snapshots_freq[plot_nodeid_x, plot_nodeid_y, :]), linewidth=2)
plot_nodeid_x = 10
plot_nodeid_y = 10
plt.semilogy(freq, np.abs(
    snapshots_freq[plot_nodeid_x, plot_nodeid_y, :]), linewidth=2)
plt.xlabel('Frequency (in Hz)')
plt.ylabel('Plate transverse deformation (in m)')
plt.title(f'FRF at node id ({plot_nodeid_x},{plot_nodeid_y}) Hz')

datadims = snapshots_freq.shape
print('dimension of the data: ', datadims)
mtrzddata = snapshots_freq.reshape((-1, datadims[-1]))
if not np.allclose(mtrzddata[:, 0], snapshots_freq[:, :, 0].reshape((-1, ))):
    raise UserWarning('this is not the right matricization')

Uk = gpu.get_podmats(sol=mtrzddata, poddim=55, plotsvs=True)
Uks = gpu.get_podmats(sol=mtrzddata, poddim=25, plotsvs=True)

rddata = Uks @ (Uks.T @ mtrzddata)
rdfrdt = rddata.reshape(datadims)

plot_nodeid_x = 0
plot_nodeid_y = 0
plt.figure()
plt.semilogy(freq, np.abs(
    rdfrdt[plot_nodeid_x, plot_nodeid_y, :]), linewidth=2)
plot_nodeid_x = -1
plot_nodeid_y = -1
plt.semilogy(freq, np.abs(
    rdfrdt[plot_nodeid_x, plot_nodeid_y, :]), linewidth=2)
plot_nodeid_x = 10
plot_nodeid_y = 10
plt.semilogy(freq, np.abs(
    rdfrdt[plot_nodeid_x, plot_nodeid_y, :]), linewidth=2)
plt.xlabel('Frequency (in Hz)')
plt.ylabel('Plate transverse deformation (in m)')
plt.title('Reduced (k=25) FRF')

# plot the deformed surface at a frequency point
plot_freqstep = 0
plt.figure()
plt.imshow(rdfrdt[:, :, plot_freqstep])
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Deformation plot at frequency {freq[plot_freqstep]} Hz')

plt.figure()
plt.imshow(snapshots_freq[:, :, plot_freqstep])
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'(k=25) Deformation plot at frequency {freq[plot_freqstep]} Hz')
plt.show()
