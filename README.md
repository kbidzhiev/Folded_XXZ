# Folded XXZ

`Folded_XXZ` is a C++ tensor-network simulation of a folded spin-1/2 XXZ
chain with three-site interactions. It uses [ITensor](https://itensor.org/)
to prepare a thermally biased state through imaginary-time evolution, then
evolve it in real time with Trotter gates.

The program is currently a single executable built from `3siteHam.cc` and
`observables.cc`.

## Requirements

- A C++ compiler compatible with the installed ITensor version.
- ITensor for C++ and its Makefile configuration.
- GNU Make.

This repository expects the local ITensor build fragments `this_dir.mk` and
`options.mk`. They are intentionally ignored by Git, so a fresh checkout will
not build until they are supplied from an ITensor installation or local setup.

## Build

After configuring the ITensor Makefile fragments expected by `Makefile`, run:

```bash
make
```

The resulting executable is `3siteHam.exe`.

## Run

Parameters are supplied as whitespace-separated `name value` pairs. All
values are numeric; unknown parameter names terminate the program after it
prints the supported parameter set.

```bash
./3siteHam.exe N 40 T 20 tau 0.01 TL 100 TR 5 EnergyProf 0.1 Sz 0.1
```

`N` is the number of sites in the original physical chain. Internally the
program uses `2 * N` spin sites: odd indices are physical sites and even
indices are their ancillas. The initial MPS pairs each physical site with its
ancilla in a singlet-like state.

The default run has `N = 10` and `T = 0`, so it prepares the initial thermal
state but does not perform positive-time real-time evolution.

## Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `N` | `10` | Physical-chain length; the MPS has `2*N` sites. |
| `J` | `1.0` | Interaction scale. |
| `tau` | `0.01` | Real-time Trotter step. |
| `T` | `0` | Final real time. |
| `dbeta` | `0.01` | Imaginary-time (inverse-temperature) step. |
| `TL`, `TR` | `100`, `5` | Left and right temperatures used for state preparation. |
| `hL`, `hR` | `0`, `0` | Left/right staggered-field amplitudes during thermal preparation. |
| `PBC` | `0` | Use periodic boundary terms when nonzero. |
| `max_bond` | `4000` | Maximum MPS bond dimension. |
| `trunc` | `1e-10` | Real-time truncation cutoff. |
| `trunc0` | `1e-10` | Imaginary-time truncation cutoff. |
| `Entropy` | `0` | Enable center-bond entropy output when nonzero. |
| `SVD_spec` | `0` | Interval for center-bond singular-value output; requires `Entropy`. |
| `Eprof` | `0` | Interval for full entanglement-entropy profiles. |
| `Energy_beta` | `1` | Enable energy versus inverse-temperature output when positive. |
| `EnergyProf` | `0` | Interval for energy and `Q1minus` profiles after preparation. |
| `Q2Prof` | `0` | Enable `Q2` profile output after preparation. |
| `Sz` | `0` | Interval for magnetization profiles after preparation. |

Several legacy parameters are accepted but are not used by the current main
simulation path: `end`, `CurrentProf`, `Current`, `XXZ`, `beta`, `write_wf`,
`TrotterOrder`, `antal`, `energy`, and `sweeps`.

For time-resolved profile options, choose intervals compatible with `tau`.
The program determines whether to write a profile using integer step counts
derived from `interval / tau`.

## Outputs

Enabled observables are written to the current working directory. Existing
files with these names are overwritten.

| File | Controlled by |
| --- | --- |
| `Entropy_center.dat` | `Entropy` |
| `SVD_spec.dat` | `SVD_spec` (and `Entropy`) |
| `Energy_beta.dat` | `Energy_beta` |
| `Entropy_profile.dat` | `Eprof` |
| `Energy_profile.dat` | `EnergyProf` |
| `Q1minus_profile.dat` | `EnergyProf` |
| `Q2_profile.dat` | `Q2Prof` |
| `Sz_profile.dat` | `Sz` |
| `Sz_average_profile.dat` | `Sz` |

Profile files separate successive time slices with blank lines, making them
convenient to plot with tools such as gnuplot.

## Code map

- `3siteHam.cc`: parameters, folded three-site Hamiltonian, Trotter evolution,
  state preparation, and output.
- `observables.h` / `observables.cc`: entropy, magnetization, energy, and
  conserved-charge measurement routines.
- `profile.h`: small scoped timing utility used by `LOG_DURATION`.
- `Pictures/`: example energy-profile animations.

## Notes

The model and observable conventions are encoded directly in the source. In
particular, the Hamiltonian acts on odd (physical) sites of the folded MPS;
the even sites are ancillas used for purification. Consult the Hamiltonian
construction in `3siteHam.cc` before comparing normalization or signs with a
different XXZ convention.
