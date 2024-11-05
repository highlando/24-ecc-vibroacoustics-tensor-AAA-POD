# from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt
# import multidim_galerkin_pod.gen_pod_utils as gpu
from multidim_galerkin_pod import ten_sor_utils as tsu

freq = np.linspace(10, 400, 391)
snapshot_data = np.load('../data/01_data.npz')
snapshots_freq = snapshot_data['snapshots_freq']

datadims = snapshots_freq.shape
print('dimension of the data: ', datadims)

dttnsr = snapshots_freq

spacedims = [2*x for x in range(10, 0, -1)]
freqqdims = [2*x for x in range(1, 30)]

spaceerr = []
spcesize = []

for dxy in spacedims:
    frqerr = []
    frqsiz = []
    cmpfrq = True
    for df in freqqdims:
        frqsiz.append(dxy*dxy*df + 2*dxy*21 + df*391)
        if df > dxy**2:
            df = dxy**2
            if not cmpfrq:
                frqerr.append(reler)
                continue
            else:
                cmpfrq = False
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
        reler = np.linalg.norm(difftens)
        print(f'crtns: {dxy:02d}x{dxy:02d}x{df:03d}: relerr = {reler:.4e}')
        frqerr.append(reler)
    spaceerr.append(frqerr)
    spcesize.append(frqsiz)
extent = freqqdims[0], freqqdims[-1], spacedims[-1], spacedims[0]
errmat = np.array(spaceerr)
plt.figure(figsize=(10, 3))
plt.imshow(np.log10(errmat), extent=extent)
plt.title('Logarithm of approximation error')
plt.xlabel('$d_f$ -- modes in frequency dimension')
plt.ylabel('$d_x = d_y$ -- modes in spatial dimension')
plt.colorbar()
plt.savefig('relapproxerr.png')
plt.figure(figsize=(10, 3))
plt.imshow(np.log10(np.array(spcesize)), extent=extent, cmap=plt.cm.cividis)
plt.colorbar()
plt.title('Logarithm of data size')
plt.xlabel('$d_f$ -- modes in frequency dimension')
plt.ylabel('$d_x = d_y$ -- modes in spatial dimension')
plt.savefig('datasize.png')
plt.show()
