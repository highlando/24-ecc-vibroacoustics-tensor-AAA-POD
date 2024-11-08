from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt
# import multidim_galerkin_pod.gen_pod_utils as gpu
from multidim_galerkin_pod import ten_sor_utils as tsu

spacepodonly = False

freq = np.linspace(10, 400, 391)
snapshot_data = np.load('../data/01_data.npz')
snapshots_freq = snapshot_data['snapshots_freq']

datadims = snapshots_freq.shape
print('dimension of the data: ', datadims)

dttnsr = snapshots_freq
nrm_snapshots_freq = np.linalg.norm(snapshots_freq)

difftl = []
lbllst = []
datadict = {'full': dttnsr}
for _fac in [2, 4, 6]:
    pdone = _fac*3
    pdtwo = _fac*3
    pdthr = _fac*_fac*3

    svdvecs_o, dttnsr_o = tsu.modek_svd(dttnsr, svddim=pdone, mode=1,
                                        return_reduced_tensor=True,
                                        plot_svs=False)
    svdvecs_ot, dttnsr_ot = tsu.modek_svd(dttnsr_o, svddim=pdtwo, mode=2,
                                          return_reduced_tensor=True,
                                          plot_svs=False)
    if not spacepodonly:
        (svdvecs_ott,
         dttnsr_ott) = tsu.modek_svd(dttnsr_ot, svddim=pdthr, mode=3,
                                     return_reduced_tensor=True, plot_svs=True)

        prjtnsr_ot = tsu.inflate_modek(dttnsr_ott, ksvecs=svdvecs_ott, mode=3)
    else:
        prjtnsr_ot = dttnsr_ot

    prjtnsr_o = tsu.inflate_modek(prjtnsr_ot, ksvecs=svdvecs_ot, mode=2)
    prjtnsr_ = tsu.inflate_modek(prjtnsr_o, ksvecs=svdvecs_o, mode=1)

    difftens = (snapshots_freq - prjtnsr_)
    difftl.append(difftens)
    if not spacepodonly:
        lbllst.append(f'{pdone}x{pdtwo}x{pdthr}')
    else:
        podstr = f'pod{pdone}x{pdtwo}'
        lbllst.append(f'{pdone}x{pdtwo}')
        datadict.update({podstr: {'dttnsr_ot': dttnsr_ot,
                                  'svdvecs_o': svdvecs_o,
                                  'svdvecs_ot': svdvecs_ot}})

(nx, ny, _) = snapshots_freq.shape

plt.figure(330, figsize=(6, 3))
clrs = [.5]
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.plasma(clrs))

# for idx in range(0, nx, 2):
#     for idy in range(0, ny, 2):
#         plt.semilogy(freq, np.abs(snapshots_freq[idx, idy, :]),
#                      linewidth=2, alpha=.1)
plt.semilogy(freq, np.linalg.norm(snapshots_freq, axis=(0, 1)),
             linewidth=2)
plt.xlabel('Frequency (in Hz)')
plt.ylabel('Plate transverse deformation (in m)')
plt.title('norm of FRF over the frequency range')
plt.tight_layout()
plt.savefig('FRF.png')

plt.figure(331, figsize=(6, 3))
clrs = [.9, .5, .2]
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(clrs))

plt.figure(331, figsize=(6, 3))

for kkk, difftens in enumerate(difftl):
    lbl = lbllst[kkk]
    plt.semilogy(freq, np.linalg.norm(difftens, axis=(0, 1)),
                 linewidth=2, alpha=.8, label=lbl)

plt.legend()
plt.xlabel('Frequency (in Hz)')
plt.ylabel('Plate transverse deformation (in m)')
plt.title('Approximation error by HOSVD of size ' +
          '$d_x \\times d_y \\times d_f$')
plt.tight_layout()
plt.savefig('FRF-reldiff-HoSVDselection.png')
savemat('podtnsrdata.mat', datadict)

plt.show()
