Vibroacoustic Data Compression and Interpolation
---

## HOSVD

 * `python/poddimcheck.py` -- checks the approximation error vs. the tensor
   reduced tensor dimensions and plots the corresponding data sizes (=number of
   entries in the core tensor and the basis vectors)
 * `python/pod_3D.py` -- computes a few HOSVD approximations

## AAA

 * `matlab/test_2D_plate_block_AAA.m` -- computes the AAA interpolation for
   various sizes of the (spatially reduced) tensor, inflates by multiplying with
   the basis vectors, and computes the difference to the full tensor.

## Data

 * `data/01_data.npz` -- full order tensor data from the simulation
 * `data/podtnsrdata.mat` -- some reduced tensors for the AAA interpolation

## Not Included

 * the python model
