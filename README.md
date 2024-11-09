Vibroacoustic Data Compression and Interpolation
===

This is the code that accompanies the paper

> Adaptive Rational Interpolation and Higher-order SVD for Low-rank Tensor Approximation in Vibroacoustics Simulation

which we (Jan Heiland, Victor Gosea, Ulrich R&ouml;mer, Davide Pradovera, Harikrishnan Sreekumar, Sabine Langer) submitted for presentation at the ECC-2025.

## Reproduction of results

To reproduce the results figures proceed as follows.

### Figure 1

This code is not provided. Please contact Ulrich R&ouml;mer (TU Braunschweig).

### Figure 2 / 3

```sh
cd python
python3 pod_3D.py  # plots of Figure 2
python3 poddimcheck.py  # plots of Figure 3
```

### Figure 4

In `matlab` run

```
cd matlab
addpath('nmodeproduct')
AAA_tensor_approximation_omega
```

### Table II/III

1. In `matlab` run

```
cd matlab
clsp=2; AAA_tensor_approximation_omega  % the values for nx=ny=6
clsp=3; AAA_tensor_approximation_omega  % the values for nx=ny=12
```

This gives the AAA errors and orders as in **Table II**.

2. Then use the script in `misc/blockdegree2hosvd.py` to calculate the data sizes
and corresponding reduced (HOSVD) frequency dimensions `tnf`. This gives that
data of **Table III**

3. Finally, run 

```sh
cd python
python3 podAAAcheck.py
```

with the `tnf`s provided, to get the *HOSVD* errors of **Table II**.

## Data

 * `data/01_data.npz` -- full order tensor data from the simulation
 * `data/podtnsrdata.mat` -- the reduced tensors for the AAA interpolation

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
