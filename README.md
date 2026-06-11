# Introduction
This is the README file for "A computational model of the mammalian auditory periphery with a closed-loop medial olivocochlear reflex simulating across-channel efferent gain control" by Guest, Farhadi, and Carney (2026), published in JASA.
This code repository focuses on replicating the figures and analyses from the paper; if you wish to use the model described in the article, instead navigate to <https://github.com/guestdaniel/EfferentModel_GuestFarhadiCarney2026>.

Please direct questions and inquiries to <daniel_guest@urmc.rochester.edu>.

# Organization
## Repository
This repository is version controlled by Git and is maintained at <https://github.com/guestdaniel/GuestFarhadiCarney2026>.
This repository is a fresh product created alongside the paper release based on older code tracked in private repositories.
As such, it does not have a complete Git version history.
Please contact the authors with questions about change histories.

## Files
```
_pkg_signal            # Local copy of AuditorySignalUtils.jl package
_pkg_utils             # Local copy of Utilities.jl internal package
_pkg_model             # Local copy of model code package (wrappers and source)
figs/                  # destination for fig files (svg, png)
src/
  ├── WarrenLiberman/  # folder of code for Warren '89 data handling
  ├── figures.jl       # figure-specific code
  ├── PrecursorRLF.jl  # custom code for RLFs with precursors
  ├── run.jl           # main script to reproduce paper
  ├── simulations.jl   # supporting code for simulations
  ├── src.jl           # package source file
  ├── stimuli.jl       # stimulus generation code
  └── supp_figures.jl  # supplemental figure code
REAMDE.md              # This markdown file
```

Note that this code depends on several mostly internal packages.
We provide frozen copies of those packages at the state they were in when this code was released, so that this code can be used standalone without worrying about other custom dependencies.
Users should not have to worry about these packages, as the main package (GuestFarhadiCarney2026) handles these issues for you (see Installation below).

# Installation and usage (paper code)
The paper code (i.e., code needed to reproduce the analyses and figures reported in the paper) is all managed by this Julia package (GuestFarhadiCarney2026).
Follow the instructions below to run the paper code (takes several hours to run).
- Open the Julia REPL in this folder.
- Press `[` in the REPL to activate the Pkg REPL.
- Type `activate .` to activate this environment.
- Type `instantiate` in the Pkg REPL to set up the environment.
- (Possibly, set up the model code, using the instructions below).
- Run `run.jl` in the Julia REPL.

# Installation and usage (model code)
The model code consists of two parts.
First, the underlying model is implemented in several C files (model source).
Second, the model can be accessed in higher-level languages (MATLAB, Julia) via convenient functions (model wrapper).

Currently, the auditory efferent model described in the paper requires a little bit of manual work to get running.
The instructions below should suffice for running this code to replicate figures and analyses from the paper.
If you wish to use the model in your own work or explore its source code, instead navigate to <https://github.com/guestdaniel/EfferentModel_GuestFarhadiCarney2026>, where the up-to-date version of the model is maintained.

Please follow the instructions below and reach out via email with any issues you encounter.

## C (model source)
The local copy of the model source code lives at `./_pkg_model/src/model`.
The permanent version is version controlled online in two repositories, <https://github.com/guestdaniel/EfferentModelCarneyC> and <https://github.com/guestdaniel/EfferentModel_GuestFarhadiCarney2026>.

You will need to compile the model for your platform.
The script `./_pkg_model/src/model/compile.sh` demonstrates how to do this using the GCC compiler.
Your tooling may vary slightly, but the goal is to compile `complex.c`, `sfie.c`, `adaptation.c`, and `model.c` into the single shared object library (`libgfc2023.so`).
Unlike the earlier variants of the (afferent) model, this model code has not seen many other users yet, and so it is expected that there will be some kinks to work out with compiling on different platforms — please reach out for assistance.

## Julia (model wrapper)
Access to the model code is provided in Julia via a model wrapper package entitled "Helios" (placeholder name).
Helios is the package contained locally at `./_pkg_model`
Follow the instructions contained therein to install and use the Helios package to run model responses.
Once installed, invoking the package with `using Helios` in a Julia script will give you access to model wrapper functions like `sim_gfc2023`. 

## MATLAB (model wrapper)
Access to the model code is provided in MATLAB via a model wrapper implemented in Mex.
The local copy of the associated code is located at `./_pkg_model/src/mex`. 
Follow the instructions contained therein to compile and use the Mex wrapper.
Once compiled and added to your MATLAB path, you will have access to model wrapper functions like `sim_efferent_model`. 

# Acknowledgements
This research was funded by NIH NIDCD R01-DC010813 (L.H.C.), NIH NIDCD F32-DC022143 (D.R.G.), and NIH NIDCD F32-DC022782 (A.F.).