import numpy as np
# requires `multidim_galerkin_pod==1.1.2`
from multidim_galerkin_pod import ten_sor_utils as tsu
# requires `multidim_galerkin_pod==1.1.2`

from helper_utils import Timer

data = np.load('./02_data_20p_test_Elognormal_5d.npz')
snapshots = data['snapshots']
print(snapshots.shape)

spcparuncdims = snapshots.shape
(xdim, ydim, pdim) = spcparuncdims[:3]
uncdims = spcparuncdims[3:]

testmode = 5
testrd = 4

psvecs, rdktnsr = tsu.modek_svd(snapshots, svddim=testrd, mode=testmode,
                                return_reduced_tensor=True)

prjtnsr = tsu.inflate_modek(rdktnsr, ksvecs=psvecs, mode=testmode)

prjer = np.linalg.norm(prjtnsr-snapshots)
fullval = np.linalg.norm(snapshots)

print(f'mode-k={testmode}: rd-k={testrd}/{spcparuncdims[testmode-1]}:' +
      f' rel-error={prjer/fullval:.4e}')

pdone = 12
pdtwo = 32

# ## Full-SVD
svddim = pdone*pdtwo

svdmode = 4
fsvdX, fsvdims = tsu.flatten_mode(snapshots, mode=4,
                                  tdims=spcparuncdims, howmany=4)

with Timer(name=f'full svd-{svddim}'):
    svdvecs, svdrtnsr = tsu.modek_svd(fsvdX, svddim=svddim, mode=svdmode,
                                      return_reduced_tensor=True)

prjtnsr = tsu.inflate_modek(svdrtnsr, ksvecs=svdvecs, mode=svdmode)

prjer = np.linalg.norm(prjtnsr-fsvdX)
fsvdatasize = svdrtnsr.size + svdvecs.size

print(f'full svd-{svddim}/{fsvdX.shape[-1]}' +
      f' rel-error={prjer/fullval:.4e}')
print(f'full svd-{svddim}/{fsvdX.shape[-1]}' +
      f' datasize={fsvdatasize}')


# ## Partial SVD

svdmodeone = 7
svdmodetwo = 4

with Timer(name=f'partial svd-{pdtwo}-{pdone}'):
    fsvdXo, fsvdso = tsu.flatten_mode(snapshots, mode=svdmodeone,
                                      tdims=spcparuncdims, howmany=1)
    fsvdXot, fsvdsot = tsu.flatten_mode(fsvdXo, mode=svdmodetwo,
                                        tdims=fsvdso, howmany=2)

    svdvecso, svdrtnsro = tsu.modek_svd(fsvdXot, svddim=pdone, mode=5,
                                        return_reduced_tensor=True)

    svdvecst, svdrtnsrot = tsu.modek_svd(svdrtnsro, svddim=pdtwo, mode=4,
                                         return_reduced_tensor=True)

prjtnsrt = tsu.inflate_modek(svdrtnsrot, ksvecs=svdvecst, mode=4)
prjtnsrot = tsu.inflate_modek(prjtnsrt, ksvecs=svdvecso, mode=5)

prjer = np.linalg.norm(prjtnsrot-fsvdXot)

print(f'partial svd-{pdtwo}/{fsvdXot.shape[-2]}-{pdone}/{fsvdXot.shape[-1]}' +
      f' rel-error={prjer/fullval:.4e}')

psvddatasize = svdrtnsrot.size + svdvecso.size + svdvecst.size

print(f'partial svd-{pdtwo}/{fsvdXot.shape[-2]}-{pdone}/{fsvdXot.shape[-1]}' +
      f' datasize={psvddatasize}')
