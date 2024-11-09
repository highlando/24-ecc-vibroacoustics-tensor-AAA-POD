# from scipy.io import savemat
import numpy as np
# import multidim_galerkin_pod.gen_pod_utils as gpu
from multidim_galerkin_pod import ten_sor_utils as tsu

freq = np.linspace(10, 400, 391)
snapshot_data = np.load('../data/01_data.npz')
snapshots_freq = snapshot_data['snapshots_freq']

datadims = snapshots_freq.shape
print('dimension of the data: ', datadims)

dttnsr = snapshots_freq

spacedimsl = [6, 12]
freqqdimsl = [[4, 5, 5, 7], [12, 14, 15, 16]]

spaceerr = []
spcesize = []

for (dxy, freqqdims) in zip(spacedimsl, freqqdimsl):
    for df in freqqdims:
        pdone = dxy
        pdtwo = dxy
        pdthr = df

        svdvecs_o, dttnsr_o = tsu.modek_svd(dttnsr, svddim=pdone, mode=1,
                                            return_reduced_tensor=True,
                                            plot_svs=False)
        svdvecs_ot, dttnsr_ot = tsu.modek_svd(dttnsr_o, svddim=pdtwo, mode=2,
                                              return_reduced_tensor=True,
                                              plot_svs=False)
        (svdvecs_ott,
         dttnsr_ott) = tsu.modek_svd(dttnsr_ot, svddim=pdthr, mode=3,
                                     return_reduced_tensor=True,
                                     plot_svs=False)

        prjtnsr_ot = tsu.inflate_modek(dttnsr_ott, ksvecs=svdvecs_ott, mode=3)

        prjtnsr_o = tsu.inflate_modek(prjtnsr_ot, ksvecs=svdvecs_ot, mode=2)
        prjtnsr_ = tsu.inflate_modek(prjtnsr_o, ksvecs=svdvecs_o, mode=1)

        difftens = (snapshots_freq - prjtnsr_)
        aprxer = np.linalg.norm(difftens)
        print(f'crtns: {dxy:02d}x{dxy:02d}x{df:03d}: err = {aprxer:.2e}')
