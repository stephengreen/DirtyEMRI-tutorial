#!/usr/bin/env python
"""Build the vendored FastEMRIWaveforms fork (the custom ``KerrCircFlux`` ODE) in
place, then make ``import few`` work from anywhere (e.g. the notebook at the repo
root).

One script, two supported environments:

* **conda / pixi** — if ``$CONDA_PREFIX`` is set, GSL/LAPACK/HDF5 are taken from
  that env and an rpath is baked in so the compiled extensions find the libs at
  runtime.
* **system libs (e.g. Google Colab)** — if there is no conda env, the build uses
  the system GSL/LAPACK/HDF5 (install them with apt first; see the Colab cell /
  README). Ubuntu's ``libgsl-dev`` is 2.7.x, which is exactly what this fork needs.

The build is CPU-only (no CUDA) and single-threaded (``--no_omp``); the trajectory
runs in well under a second, so OpenMP buys nothing here and ``--no_omp`` keeps the
build uniform across macOS (Apple clang rejects bare ``-fopenmp``) and Linux.
"""
import os
import sys
import sysconfig
import pathlib
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
FEW = HERE / "FastEMRIWaveforms"

if not (FEW / "setup.py").exists():
    sys.exit("ERROR: cannot find FastEMRIWaveforms/setup.py next to this script.")

env = os.environ.copy()
# setup.py does `import distutils`; on Python >=3.12 that only resolves via
# setuptools' vendored copy. Harmless on older Pythons.
env.setdefault("SETUPTOOLS_USE_DISTUTILS", "local")

cmd = [sys.executable, "setup.py", "build_ext", "--inplace", "--no_omp"]

prefix = os.environ.get("CONDA_PREFIX")
if prefix:
    prefix = pathlib.Path(prefix)
    libdir, incdir = prefix / "lib", prefix / "include"
    cmd += ["--gsl", str(prefix), "--lapack", str(prefix)]
    env["LDFLAGS"] = f"-L{libdir} -Wl,-rpath,{libdir} " + env.get("LDFLAGS", "")
    env["CPATH"] = str(incdir) + (os.pathsep + env["CPATH"] if env.get("CPATH") else "")
    print(f"[build_few] using conda/pixi env at {prefix}")
else:
    print("[build_few] no $CONDA_PREFIX -> using system GSL/LAPACK/HDF5 "
          "(make sure libgsl-dev, liblapacke-dev, libhdf5-dev are installed)")

print("[build_few] >>", " ".join(cmd))
subprocess.check_call(cmd, cwd=str(FEW), env=env)

# Register the build so `import few` AND `from pyInspiral import ...` resolve from
# any working directory. The compiled extensions are top-level modules that land
# in FastEMRIWaveforms/, so we put that directory on the path via a .pth file in
# the active environment's site-packages.
purelib = pathlib.Path(sysconfig.get_paths()["purelib"])
pth = purelib / "few_tutorial.pth"
try:
    pth.write_text(str(FEW) + "\n")
    print(f"[build_few] wrote {pth} -> {FEW}")
except OSError as e:
    print(f"[build_few] WARNING: could not write {pth} ({e}); "
          f"the notebook will add the path itself.")

# Smoke test in a fresh interpreter (picks up the .pth automatically).
smoke = (
    "from few.trajectory.inspiral import EMRIInspiral;"
    "traj=EMRIInspiral(func='KerrCircFlux');"
    "r=traj(1e6,50.0,0.9,15.482608237080893,0.0,1.0,1e-5,8.0,0.0,T=1.0,dt=10.0);"
    "print('[build_few] SMOKE TEST OK: KerrCircFlux trajectory has',len(r[0]),'points')"
)
subprocess.check_call([sys.executable, "-c", smoke], cwd=str(HERE), env=env)
print("[build_few] FastEMRIWaveforms (KerrCircFlux) built successfully.")
