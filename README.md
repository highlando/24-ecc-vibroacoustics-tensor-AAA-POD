Vibroacoustic Data Compression and Interpolation
---

This is the code that accompanies the paper

> Adaptive Rational Interpolation and Higher-order SVD for Low-rank Tensor Approximation in Vibroacoustics Simulation

which we (Jan Heiland, Victor Gosea, Ulrich R&ouml;mer, Davide Pradovera, Harikrishnan Sreekumar, Langer Sabine) submitted for presentation at the ECC-2025.

To reproduce the results figures proceed as follows.

### Figure 1

This code is not provided. Please contact Ulrich R&ouml;mer (TU Braunschweig).

### Figure 2 / 3

```sh
cd python
python3 pod_3D.py  # plots of Figure 2
python3 poddimcheck.py  # plots of Figure 3
```

## Figure 4

In `matlab` run

```
cd matlab
addpath('nmodeproduct')
AAA_tensor_approximation_omega
```



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

## Dependencies

### Python

Use/see the `python/requirements.txt`

```sh
pip install -r requirements.txt`
```

Among others, we need

 * `multidim_galerkin_pod` -- required >= 1.1.3 / tested with 1.1.3
 * `matplotlib` -- tested with 3.9.2
 * `numpy` -- tested with 2.1.3
 * `scipy` -- tested with 1.14.1
