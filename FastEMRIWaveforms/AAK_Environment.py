import time
import numpy as np
import warnings
import matplotlib.pyplot as plt
try:
    import cupy as xp

    gpu_available = True

except (ModuleNotFoundError, ImportError) as e:
    import numpy as xp

    warnings.warn(
        "CuPy is not installed or a gpu is not available. If trying to run on a gpu, please install CuPy."
    )
    gpu_available = False


from few.trajectory.inspiral import EMRIInspiral
from few.utils.baseclasses import SchwarzschildEccentric, Pn5AAK, ParallelModuleBase
from few.waveform import AAKWaveformBase, Pn5AAKWaveform
from few.summation.aakwave import AAKSummation
from few.utils.utility import get_mismatch

insp_kwargs = {
    "err": 1e-10,
    "DENSE_STEPPING": 0,
    "max_init_len": int(1e4),
    "func":"KerrCircFlux"
    }

# keyword arguments for summation generator (AAKSummation)
sum_kwargs = {
    "use_gpu": gpu_available,  # GPU is availabel for this type of summation
    "pad_output": False,
}

num_threads = 1

waveform_class = AAKWaveformBase(
            EMRIInspiral,
            AAKSummation,
            inspiral_kwargs=insp_kwargs,
            sum_kwargs=sum_kwargs,
            use_gpu=gpu_available,
            num_threads=num_threads,
        )

wave_generator = Pn5AAKWaveform( sum_kwargs=sum_kwargs, use_gpu=gpu_available)

####################
# set initial parameters
M = 1e6
mu = 50.0
p0 = 15.482608237080893
e0 = 0.0
a = 0.9
x0 = 1.0
Y0 = 1.0
# accretion parameters from https://arxiv.org/abs/2207.10086
A = 1e-5
n = 8.0
m = 0.0

qS = 0.2
phiS = 0.2
qK = 0.8
phiK = 0.8
dist = 1.0
mich = False
dt = 10.0
T = 1.0

if gpu_available:
    tic = time.perf_counter()
    for _ in range(10):
        output1 = waveform_class(M,mu,a,p0,e0,Y0,dist,qS,phiS,qK,phiK,A,n,m,
            Phi_phi0=0.0,Phi_theta0=0.0,Phi_r0=0.0,mich=False,dt=dt,T=T)
    toc = time.perf_counter()
    print((toc-tic)/10)

# check limit
output1 = waveform_class(M,mu,a,p0,e0,Y0,dist,qS,phiS,qK,phiK,A,n,m,
        Phi_phi0=0.0,Phi_theta0=0.0,Phi_r0=0.0,mich=False,dt=dt,T=T).get()

A = 0.0
output2 = waveform_class(M,mu,a,p0,e0,Y0,dist,qS,phiS,qK,phiK,A,n,m,
        Phi_phi0=0.0,Phi_theta0=0.0,Phi_r0=0.0,mich=False,dt=dt,T=T).get()
print(output2)
print(get_mismatch(output2.real, output1.real))

plt.figure()
plt.plot(output1.real[-1000:],label='A non zero')
plt.plot(output2.real[-1000:],'--',label='A = 0')
plt.legend()
plt.show()
