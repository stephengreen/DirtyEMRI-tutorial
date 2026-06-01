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
