### Rocket Engine Simulation 🚀

This project is a continuation of Ben Klammer's Mech 498 Thesis on Hybrid Modelling [1]. 

This program was designed to be a modular framework for a user to combine and use different models to simulate rocket engines. This decision to make the program modular will hopefully allow us to add more detail if/as required throughout the engine design cycle.

The modelling is especially focused on nitrous oxide as a self-pressurant.


## RUNNING THE PROGRAM:
This file takes input in the form of a python module file. 

Ensure the terminal is in the directory this file is in:  #python3 -m src <filename of a desired python module file without the .py>

The inputs folder contains the input files where the user can update the inputs the program uses.

The analysis_mode list in the input file selects which models are active.

The user interface isn't great, but if the user interface is good, are you really innovating?

### Running Paraffin Hybrids --> See the instructions in PARAFFIN_DEFINITION.md

## Getting Started (instructions for windows)

This project is intended to run inside WSL (Windows Subsystem for Linux) with a Python virtual environment (venv). This is because I was only able to get the rocketCEA package to work here. Allegedly there have been some changes that make this install easier since then, but I have not tried them.

Setup Steps:

Install WSL

Microsoft Guide: https://learn.microsoft.com/en-us/windows/wsl/install

Recommended distro: Ubuntu 22.04

Install the fortran compiler (Ubuntu instructions): https://rocketcea.readthedocs.io/en/latest/installgfortran.html

Install Python (via pyenv)

Guide: https://github.com/pyenv/pyenv#installation

Tested with Python 3.11.x

Clone the repo
>git clone <repo_url>
>cd <repo>

Create and activate a virtual environment
python3 -m venv venv
source venv/bin/activate

Install required libraries
A requirements.txt file is included so everyone uses the same library versions.
>pip install -r requirements.txt

To see what’s currently installed in your venv:
>pip list

To update requirements.txt after adding a new library:
>pip freeze > requirements.txt

Always create a new branch for edits (git checkout -b <branch_name>).


## Sources and Citations:
| Number | Source                            | Contribution / Use Description                 | Link to Source    |
|--------|-----------------------------------|------------------------------------------------|-------------------|
| [1]    | NASA Chemical Equilibrium with Applications (CEA) | Used with the rocketcea library for combustion calculations | https://cearun.grc.nasa.gov/ |
| [2]    | Ben Klammer, Hybrid Modelling: https://github.com/bklammer/HybridModeling | Translated equilibrium model used in this thesis into python and used it for analysis |https://github.com/bklammer/HybridModeling |
| [3]    | Benjamin S. Waxman, Jonah E. Zimmerman, Brian J. Cantwell, Mass Flow Rate and Isolation Characteristics of Injectors for Use with Self-Pressurizing Oxidizers in Hybrid Rockets | Used to implement HEM and Dyer injector models, experimental data to validate Emerson + Mohammad injector model and understand nitrous physics. |https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf |
| [4]    | Numerical Modeling of Pressurization of a Propellant Tank | Referenced for setting up the pressurized fuel tank model, although my model is much more simple at this stage | https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-tankpressurization-jpp2001.pdf?emrc=66201987b6c8c |
| [5]    | Emerson Vargas Niño, Mohammad Reza H. Razavi, Design of Two-Phase Injectors Using Analytical and Numerical Methods with Application to Hybrid Rockets | Used to implement Emerson + Mohammad injector model | https://emersonvn.com/project/two_phase_injector/# |
| [6]    | Tomasz Palacz, Jacek Cieslik. Experimental Study on the Mass Flow Rate of the Self Pressurizing Propellants in the Rocket Injector | Used experimental data to build a better understanding of model feed system | https://www.researchgate.net/publication/355773008_Experimental_Study_on_the_Mass_Flow_Rate_of_the_Self-Pressurizing_Propellants_in_the_Rocket_Injector  ~ DOI:10.3390/aerospace8110317 |
| [7]    | Zilliac, Gregory & Karabeyoglu, Arif. (2005). Modeling of Propellant Tank Pressurization. | Used to setup non-equilibrium tank model (1/2) | https://arc.aiaa.org/doi/10.2514/6.2005-3549 |
| [8]    | Review and Evaluation of Models for Self-Pressurizing Propellant Tank Dynamics, Jonah E. Zimmerman, Benjamin S. Waxman, Brian Cantwell and Greg Zilliac | Used to setup non-equilibrium tank model (2/2) | https://arc.aiaa.org/doi/10.2514/6.2013-4045 |
| [9]    | Michigan Aeronautical Science Association Liquid Bi-Propellant Rocket Team 37 Project Technical Report for the 2018 IREC | Script Test Case | https://www.soundingrocket.org/uploads/9/0/6/4/9064598/37_project_report.pdf |
| [10]   | Technical report for Stanford SSI's rocket Project Olympus for the Spaceport America Cup Competition 2022 | Script Test Case | https://purl.stanford.edu/tx770vm3347 |
| [11]   | Two-phase flow through pressure safety valves. Experimental investigation and model prediction, Gino Boccardi et al | explains omega model (part of [5]) | https://www.sciencedirect.com/science/article/pii/S0009250905003313 |
| [12]   | Effects of Injector Design and Impingement Techniques on the Atomization of Self-Pressurizing Oxidizers, Benjamin S. Waxman, Brian J. Cantwell, Greg Zilliac | nitrous injector | https://arc.aiaa.org/doi/10.2514/6.2012-3906 |
| [13]   | Rocket Propellant and Pressurization Systems Elliot Ring | used to help better understand pressurized tank models | very hard to find- out of print and pdf is hard to come by |
| [14]   | NIST Chemical Webbook | Used to validate thermo properties from libraries | https://webbook.nist.gov/cgi/fluid.cgi?ID=C10024972&Action=Page |
| [15]  | Thorade, M., Saadat, A. (2013): Partial derivatives of thermodynamic state properties for dynamic simulation. ‐ Environmental Earth Sciences, 70, 8, 3497‐3503 | for span wagner eos partial derivs |https://gfzpublic.gfz-potsdam.de/rest/items/item_247373_5/component/file_306833/content?download=true
| [16]  | Modeling Feed System Flow Physics for Self-Pressurizing Propellants | Understanding Metastable States | https://www.researchgate.net/publication/268482381_Modeling_Feed_System_Flow_Physics_for_Self-Pressurizing_Propellants |
| [17]  | Span, Multiparameter Equations of State An Accurate Source of Thermodynamic Property Data | Understanding how to implement the Span Wagner EOS | https://www.researchgate.net/publication/40381676_Multiparameter_equations_of_state_an_accurate_source_of_thermodynamic_property_data_with_151_figures_and_tables |
| [18]  | A Reliable and Useful Method to Determine the Saturation State from Helmholtz Energy Equations of State | Another source that helps explain the iterative algo to solve saturation properties w helmholtz eos | https://www.jstage.jst.go.jp/article/jtst/3/3/3_3_442/_pdf |
| [19]  |Mass Flow Rate and Isolation Characteristics of Injectors for Use with Self-Pressurizing Oxidizers in Hybrid Rockets | Inj Design | https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf |
| [20]  | AIAA Liquid Rocket Thrust Chambers | | |
| [21]  | Huzel + Huang | | |
| [22]  | Ben Klammer 446 CC model | | |
| [23]  | Implementing CEA calculations using Cantera - Kyle Niemeyer | Starting to outgrow rocketcea, looking into using cantera for cc model because it is more compatible with differential eqns | https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html |
| [24]  | NASA SP-1311-2 | fixing paraffin definition in rocketcea. see PARAFFIN_DEFFINITION.md | https://shepherd.caltech.edu/EDL/PublicResources/sdt/refs/NASA-RP-1311-2.pdf |
| [25]  | Computational model for performance prediction of a nitrous oxide / eicosane hybrid rocket engine | hybrid cc model that integrates P_cc_dot | https://www.researchgate.net/publication/381880279_Computational_model_for_performance_prediction_of_a_nitrous_oxide_eicosane_hybrid_rocket_engine|
| [26]  | Liquid Rocket Propulsion Instability | trying to understand instability for prelim/detail design |
| [27]  | AGARD Heat Transfer in Rocket Engines | Estimating Heat Transfer in Injector | https://apps.dtic.mil/sti/trecms/pdf/AD0733362.pdf |
| [28]  | Liquid Rocket Engine Injectors | | https://ntrs.nasa.gov/api/citations/19760023196/downloads/19760023196.pdf |


Ryan Wright, From Calgary, Alberta, Canada   

<!-- 
source venv/bin/activate 
uvr working link: https://github.com/UVicRocketry/uvr_hybrid_modelling_v2

-->