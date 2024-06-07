# parallel_energy.py
# James Widdicombe
# Last Updated 16/12/2019

# Load the modules
from matplotlib import rcParams
import matplotlib.pyplot as plt
import yt
import numpy as np
import time
import matplotlib

matplotlib.use("Agg")

# Timings
start_time = time.time()

# Enable Parallelism
yt.enable_parallelism()

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "chk/Sim_000*.3d.hdf5"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = ts[0].domain_right_edge / 2.0
adjusted_right = int(center[0]) * 2 - 8
print(center)


def _detgamma(field, data):
    return data["chi"]**(-3.0)


def _rhoJ(field, data):
    return -(-data["y"] * data["s1"] + data["x"] *
             data["s2"]) * data["detgamma"]**(1.0 / 2.0)


def _rhofixed(field, data):
    return (data["rho"] * data["lapse"] - data["s1"] * data["shift1"] - data["s2"]
            * data["shift2"] - data["s3"] * data["shift3"]) * data["detgamma"]**(1.0 / 2.0)


# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):

    i.add_field("detgamma", _detgamma, units="")
    i.add_field("rhoJ", _rhoJ, units="cm")
    i.add_field("rhofixed", _rhofixed, units="")

    # Defining domain of integration
    ad = i.r[
        int(center[0]) - 16: int(center[0]) + 16,
        int(center[0]) - 16: int(center[0]) + 16,
        0: 16
    ]
    vol = ad.sum("cell_volume")
    mass = vol * ad.mean("rhofixed", weight="cell_volume")
    angmom = vol * ad.mean("rhoJ", weight="cell_volume")

    array = [i.current_time, mass, angmom]

    sto.result = array
    sto.result_id = str(i)
if yt.is_root():
    timedata = []
    rhodata = []
    angmom = []

    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        rhodata.append(L[1][1])
        angmom.append(L[1][2])
    np.savetxt("time.out", timedata)
    np.savetxt("energy.out", rhodata)
    np.savetxt("angular_momentum.out", angmom)

    plt.figure(1)
    plt.plot(timedata, rhodata)
    plt.ylabel("Energy")
    plt.xlabel("Time $[1/m]$")
    plt.grid()
    plt.savefig("Energy.png")
    plt.close()

    plt.figure(1)
    plt.plot(timedata, angmom)
    plt.ylabel("Angular Momentum")
    plt.xlabel("Time $[1/m]$")
    plt.grid()
    plt.savefig("angular_momentum.png")
    plt.close()
