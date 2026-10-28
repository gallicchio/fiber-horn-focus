# Reproduce everything.  `make envs` once, then `make sims` (~23 h on a 16-core desktop),
# then `make media`.  Raw output goes to data/ (about 15 GB, not in git).
SIM_ENV    ?= $(HOME)/micromamba/envs/fiberlaunch
RENDER_ENV ?= $(HOME)/micromamba/envs/fiberlaunch-render
PY  := $(SIM_ENV)/bin/python
PYR := $(RENDER_ENV)/bin/python
RES ?= 16
ROOT := data/main/res$(RES)

.PHONY: all envs sims main sweep convergence bpm horn media figures movies render previews test

all: sims media

envs:
	micromamba create -y -f environment.yml
	micromamba create -y -f environment-render.yml

sims: main sweep convergence bpm horn

main:
	RES=$(RES) scripts/run_main.sh

sweep:
	RES=$(RES) scripts/run_sweep.sh

convergence:
	scripts/run_convergence.sh

bpm:
	$(PY) scripts/run_bpm.py --root $(ROOT)

horn:
	scripts/run_horn.sh

media: figures movies render previews

figures:
	$(PY) scripts/fig_theory.py
	$(PY) scripts/fig_results.py --root $(ROOT)
	$(PY) scripts/fig_scales.py --root $(ROOT)
	$(PY) scripts/summarize_convergence.py | tee media/convergence.txt
	$(PY) scripts/check_steady.py --root $(ROOT) | tee media/steady_state.txt
	$(PY) scripts/horn_design.py $(foreach r,20 40 60 80 120,data/horn/design/emit_res$(r)) | tee media/horn_design.txt
	$(PY) scripts/horn_media.py --only side modes budget

movies:
	$(PY) scripts/movie_turnon.py --root $(ROOT)
	$(PY) scripts/movie_flythrough.py --root $(ROOT)
	$(PY) scripts/media_bpm.py
	$(PY) scripts/horn_media.py --only movie

render:
	$(PY) scripts/export_volumes.py --root $(ROOT)
	$(PYR) scripts/render_3d.py

previews:
	scripts/make_previews.sh

test:
	$(PY) -m pytest -q tests
