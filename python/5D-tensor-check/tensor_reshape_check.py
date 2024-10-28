import numpy as np
from multidim_galerkin_pod import ten_sor_utils as tsu

data = np.load('./02_data_20p_test_Elognormal_5d.npz')
snapshots = data['snapshots']
print(snapshots.shape)
spcparuncdims = snapshots.shape
(xdim, ydim, pdim) = spcparuncdims[:3]
uncdims = spcparuncdims[3:]

ssf, ftdims = tsu.flatten_mode(snapshots, 1, tdims=spcparuncdims)
ssff, fftdims = tsu.flatten_mode(ssf, 3, tdims=ftdims)

ss, fdims = tsu.unflatten_mode(ssff, 3, ftdims=fftdims)

sss, sfdims = tsu.unflatten_mode(ssff, 1, ftdims=fftdims)

ssss, ssfdims = tsu.unflatten_mode(sss, 4, ftdims=sfdims)

print(ssff.shape)
print(ftdims)

print(ssf.shape)
print(fftdims)

print(ss.shape)
print(fdims)

print(sss.shape)
print(sfdims)

print(ssss.shape)
print(ssfdims)

print(np.linalg.norm(ssss-snapshots))


sf, ftdims = tsu.flatten_mode(snapshots, 6, tdims=spcparuncdims, howmany=2)
sf, ftdims = tsu.flatten_mode(sf, 4, tdims=ftdims, howmany=1)
print(sf.shape)
print(ftdims)
