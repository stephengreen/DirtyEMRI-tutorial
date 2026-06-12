# Kerr Circular EMRIs in Accretion Disks — new Tutorial

The newest version of this tutorial can be found in the notebook 
[Trajectory_tutorial_environments.ipynb](Trajectory_tutorial_environments.ipynb).
The notebook uses the latest version of FastEMRIWaveforms and customizes the trajectory 
ODE to account for interactions with the environment.

## Option A — Conda environment

Before using the notebook, create the following conda environment from your terminal:

```bash
conda env create -f few_tut_new.yml
conda activate few_tut_new
python -m ipykernel install --user --name few_tut_new --display-name "Python (few_tut_new)"
jupyter lab Trajectory_tutorial_environments.ipynb
```

## Credits

The Trajectory_tutorial_environments tutorial
is based on the [Trajectory_tutorial](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms/blob/master/examples/Trajectory_tutorial.ipynb) by Christian Chapman-Bird, Lorenzo Speri and Maxime Pigou and the [Tutorial_FrequencyDomain_Waveforms](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms/blob/master/examples/Tutorial_FrequencyDomain_Waveforms.ipynb) by Michael Katz and Lorenzo Speri, with tutorial revisions and the homework exercises by **Laura Sberna**. It builds on
the [Fast EMRI Waveforms (FEW)](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms)
framework by Michael Katz, Alvin Chua, Niels Warburton, Scott Hughes, and other
FEW developers.

## License

This project is licensed under the **GNU General Public License v3** — see
[`LICENSE`](LICENSE) for the full text and [`COPYRIGHT`](COPYRIGHT) for
authorship and attribution.

The bundled `FastEMRIWaveforms/` fork was taken from FEW 1.x in 2023, when FEW
was GPLv3. (Upstream FEW [relicensed to MIT in February 2025](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms/blob/master/LICENSE),
after this fork was made, so that newer MIT license does not apply to the
snapshot here.) The fork also includes third-party GPLv3 code (Tensor Product
Interpolation, © Michael Pürrer & Jonathan Blackman), so the repository as a
whole is distributed under the GPLv3.

# Kerr Circular EMRIs in Accretion Disks — old Tutorial

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
2. Run the **first code cell**. It installs the dependencies and then **restarts
   the runtime — this is expected and normal.**
3. After the restart, choose **Runtime → Run all.**

The second cell builds the fork (~2-3 minutes the first time) and then the physics
cells just run.

**Two things during setup look alarming but are normal:** (1) a block of red
`ERROR: pip's dependency resolver...` text about NumPy — harmless, it's just
Colab's other packages; the install succeeded. (2) The runtime **restarts** after
the first cell (possibly shown as *"session crashed"*) — intentional, for a clean
NumPy. Just **Runtime → Run all** afterward.

> Why the restart? The fork needs `numpy < 2`, and downgrading NumPy inside a live
> kernel corrupts it — so we install everything, restart once for a clean NumPy,
> then build against Colab's own Python. (Colab's Ubuntu already ships GSL 2.7, the
> version the fork needs.)

## Option B — Local, with [pixi](https://pixi.sh) (recommended for offline use)

[pixi](https://pixi.sh) builds a self-contained environment (compiler, GSL,
LAPACK, Python 3.12, Jupyter) — no conda/brew/apt and no system libraries to
install yourself.

**1. Install pixi** (once — it's a single self-contained binary):

```bash
# macOS / Linux
curl -fsSL https://pixi.sh/install.sh | sh
```
```powershell
# Windows (PowerShell)
iwr -useb https://pixi.sh/install.ps1 | iex
```

…or `brew install pixi`. It installs to `~/.pixi/bin` — **open a new terminal**
afterward so it's on your `PATH`. (Details: <https://pixi.sh/latest/#installation>.)

**2. Run the tutorial:**

```bash
git clone https://github.com/stephengreen/DirtyEMRI-tutorial.git
cd DirtyEMRI-tutorial
pixi run lab        # builds the fork (once) and opens it in JupyterLab
```

On first run, pixi downloads everything into a local `.pixi/` folder — give it a
couple of minutes. To instead just check it all builds and the notebook runs
end-to-end:

```bash
pixi run verify
```

---

## Option C — Local, with your own conda / mamba

If you already use conda or mamba and would rather manage the environment
yourself:

```bash
conda create -n few_tut -c conda-forge python=3.12 "numpy<2" "cython<3" scipy \
    gsl=2.7 lapack liblapacke openblas hdf5 h5py requests tqdm matplotlib jupyterlab
conda activate few_tut
python build_few.py        # compiles the fork against this env + registers `few`
jupyter lab Playground_KerrCircularAccretion.ipynb
```

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
├── LICENSE                                   GNU General Public License v3
├── COPYRIGHT                                 authorship & attribution
└── README.md
```

## Credits

The DirtyEMRI tutorial — the physics, the custom `KerrCircFlux` accretion-disk
trajectory, and the [`Playground_KerrCircularAccretion.ipynb`](Playground_KerrCircularAccretion.ipynb)
notebook — was created by **[Lorenzo Speri](https://github.com/lorenzsp)**, with
tutorial revisions and the homework exercises by **Laura Sberna**. It builds on
the [Fast EMRI Waveforms (FEW)](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms)
framework by Michael Katz, Alvin Chua, Niels Warburton, Scott Hughes, and the
FEW developers.

This repository repackages their tutorial to run out-of-the-box (the
cross-platform build, the pixi setup, and the self-configuring Colab cells)
— Stephen Green. The science is Lorenzo's and Laura's.

## License

This project is licensed under the **GNU General Public License v3** — see
[`LICENSE`](LICENSE) for the full text and [`COPYRIGHT`](COPYRIGHT) for
authorship and attribution.

The bundled `FastEMRIWaveforms/` fork was taken from FEW 1.x in 2023, when FEW
was GPLv3. (Upstream FEW [relicensed to MIT in February 2025](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms/blob/master/LICENSE),
after this fork was made, so that newer MIT license does not apply to the
snapshot here.) The fork also includes third-party GPLv3 code (Tensor Product
Interpolation, © Michael Pürrer & Jonathan Blackman), so the repository as a
whole is distributed under the GPLv3.
