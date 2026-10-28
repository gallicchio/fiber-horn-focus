# Shared settings for the run scripts.  Source this file; override any variable first.
#   PYTHON  the simulation environment's python (Meep, MPI build)
#   NP      MPI processes (one per physical core is best; Meep is memory-bound)
#   DATA    where raw simulation output goes (several GB; not in git)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
: "${MAMBA_ROOT_PREFIX:=$HOME/micromamba}"
: "${PYTHON:=$MAMBA_ROOT_PREFIX/envs/fiberlaunch/bin/python}"
: "${MPIRUN:=$(dirname "$PYTHON")/mpirun}"
: "${NP:=16}"
: "${DATA:=$ROOT/data}"
# Single-machine runs: use shared memory, and keep UCX from trying RDMA on any
# InfiniBand/RoCE adapters (which fails without an unlimited memlock limit).
export UCX_TLS="${UCX_TLS:-self,sm,tcp}"
export OMP_NUM_THREADS=1
cd "$ROOT"
fdtd() { "$MPIRUN" -np "$NP" "$PYTHON" -m fiberlaunch.fdtd "$@"; }
