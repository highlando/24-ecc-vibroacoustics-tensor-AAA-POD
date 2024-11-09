import numpy as np
import matplotlib.pyplot as plt
import multidim_galerkin_pod.gen_pod_utils as gpu

freq = np.linspace(10, 400, 391)
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
# plot_freqstep = 0
# plt.figure()
# plt.imshow(rdfrdt[:, :, plot_freqstep])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Deformation plot at frequency {freq[plot_freqstep]} Hz')
# plt.figure()
# plt.imshow(snapshots_freq[:, :, plot_freqstep])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'(k=25) Deformation plot at frequency {freq[plot_freqstep]} Hz')
plt.show()
