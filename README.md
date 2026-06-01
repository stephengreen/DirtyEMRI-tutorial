# Kerr Circular EMRIs in Accretion Disks — Tutorial

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/stephengreen/DirtyEMRI-tutorial/blob/main/Playground_KerrCircularAccretion.ipynb)

A slimmed-down, ready-to-run version of the **DirtyEMRI** playground for modelling
Extreme Mass Ratio Inspirals (EMRIs) of compact objects on circular Kerr orbits
embedded in an accretion disk. Everything you need to run
[`Playground_KerrCircularAccretion.ipynb`](Playground_KerrCircularAccretion.ipynb)
is in this one repository.

The notebook uses a **custom `KerrCircFlux` trajectory** that lives only in the
vendored `FastEMRIWaveforms/` fork here — it is *not* in the public PyPI
`fastemriwaveforms` package, so the fork is built from source. No GPU is required.

---

## Option A — Google Colab (nothing to install)

1. Open the notebook in Colab (the **"Open in Colab"** badge above, or the one at
   the top of the notebook).
2. Run the **first code cell**. It installs `conda` and then **restarts the
   runtime — this is expected and normal.**
3. After the restart, choose **Runtime → Run all.**

Setup picks up automatically: it installs the tested dependencies (via conda, so
the right GSL/HDF5 are used) and builds the fork — ~2-3 minutes the first time,
then the physics cells just run.

> Why conda on Colab? The fork needs GSL < 2.8 and a specific NumPy/Cython, which
> conflict with Colab's stock Python. `condacolab` gives an isolated, known-good
> environment — hence the one-time restart.

## Option B — Local, with [pixi](https://pixi.sh) (recommended for offline use)

[pixi](https://pixi.sh/latest/#installation) installs a self-contained
environment (compiler, GSL, LAPACK, Python, Jupyter) — no conda/brew/apt needed.

```bash
git clone https://github.com/stephengreen/DirtyEMRI-tutorial.git
cd DirtyEMRI-tutorial
pixi run lab        # builds the fork (once) and opens the notebook in JupyterLab
```

That's it. To just verify the environment builds and the notebook runs end-to-end:

```bash
pixi run verify
```

## Option C — Local, with your own conda/mamba

```bash
conda create -n few_tut -c conda-forge python=3.9 "numpy<2" "cython<3" scipy \
    gsl=2.7 lapack liblapacke openblas hdf5 h5py requests tqdm matplotlib jupyterlab
conda activate few_tut
python build_few.py
jupyter lab Playground_KerrCircularAccretion.ipynb
```

---

## Why these specific versions? (the install gotchas)

This fork is from ~2022 and has a few hard requirements that trip up a naive
install. They are all handled for you above, but for the record:

| Requirement | Why |
|---|---|
| **GSL < 2.8** (e.g. 2.7) | GSL 2.8 rewrote the B-spline API the fork uses; it will not compile against 2.8. (Homebrew ships 2.8 → fails. Ubuntu/Colab ship 2.7 → fine.) |
| **NumPy < 2**, **Cython < 3** | The Cython/C++ was written against the NumPy 1.x C API. |
| `--no_omp` on macOS | Apple's clang rejects bare `-fopenmp`. The build is single-threaded (plenty fast). |
| `lapacke.h` from `openblas` | On conda-forge the header ships with the `openblas` package, not `liblapacke`. |

No internet/Zenodo download is needed at runtime — the `KerrCircFlux` model reads
two small data tables already included in `FastEMRIWaveforms/few/files/`.

## Forking

The Colab badge and the notebook's `REPO_URL` point at
`github.com/stephengreen/DirtyEMRI-tutorial`. If you fork or rename the repo,
update those two spots (the badge at the top of the notebook / this README, and
`REPO_URL` in the notebook's first cell) to your own URL.

## What's inside

```
DirtyEMRI-tutorial/
├── Playground_KerrCircularAccretion.ipynb   the tutorial
├── FastEMRIWaveforms/                        the custom FEW fork (KerrCircFlux)
├── build_few.py                              cross-platform build (conda/pixi or Colab)
├── pixi.toml                                 local one-command setup
└── README.md
```
