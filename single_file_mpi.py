import itertools

import numpy as np
from mpi4py import MPI

import pyasdf

# Use default communicator. No need to complicate things.
COMM = MPI.COMM_WORLD


# File to read from - force MPI off within ASDF - we will manually do it
# outside of it. Also open in read-only mode to not accidentially mess with the
# file.
ds = pyasdf.ASDFDataSet("./asdf_example.h5", mpi=False, mode="r")


# Collect all available traces.
if COMM.rank == 0:
    # Get ALL traces. Looks awkward but its essentially a double loop.
    jobs = list(itertools.chain.from_iterable(
        (_i for _i in ds.waveforms[station].list() if _i != "StationXML")
        for station in ds.waveforms.list()))
    # Split into however many cores are available.
    jobs = [jobs[_i::COMM.size] for _i in range(COMM.size)]
else:
    jobs = None

# Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)

results = []

# This will essentially be a loop over traces.
for job in jobs:
    # Read the trace ... again slightly awkward.
    tr = ds.waveforms['.'.join(job.split('.')[:2])][job][0]

    # Now just do your thing for each trace.
    results.append(
        {"trace_name": job,
         "output": np.ones((9, 10)) * tr.data.mean()})


# Gather results on rank 0.
results = MPI.COMM_WORLD.gather(results, root=0)


if COMM.rank == 0:
    # Flatten list of lists.
    results = [_i for temp in results for _i in temp]

    # At this point you could loop over the results and write a new ASDF file
    # and write the results to the auxiliary data. Remember to set `mpi=off`
    # when initializing the ASDFDataSet().

    print("Results:", results)
