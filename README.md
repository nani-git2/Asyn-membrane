Model to study surface condensation of 𝛼-synuclein on a lipid bilayer. The membrane is composed of two kinds of lipids (DOPC and DOPS) and is represented via a 2D Ising model (lipid-lipid interaction Jm). 𝛼-synuclein proteins are represented using a 3D lattice-gas model (protein-protein interaction Jp). The proteins interact with the membrane via tethers (protein-tether interaction Jt)


The custom code is written in Julia language (Version 1.10.9) and was exectued on a standard desktop (32 GB RAM, 12th Gen Intel(R) Core(TM) i7-12700) with Debian 12.11 operating system. It was also tested on a Lenovo latop with 32BG RAM and Intel(R) Core(TM) Ultra 9 185H processor, with Ubuntu 22.04.5 LTS operating system.

To run the code, just save the .jl file to a machine with Julia installed (https://julialang.org/install/). Specify the parameters and either run 
>> output = main()
in the REPL, or include a line to save the output data in a preferred format and analyse later. For Julia code, JLD is a recommended format.

The output gives the state of the system: particle and tether positions, and the membrane states. The excess density is computed from this as mentioned in the manuscript.


