### Note: mostly broken code! do not use right now lol --> see DEV LRE Progress Update... for the most recent documentation and explaination for where this program is rn

### Rocket Engine Simulation 🚀

This program was designed to be a framework for a user to combine and use different models to simulate a rocket engine.

## RUNNING THE PROGRAM:
This file takes input in the form of a python module file. 

In the directory this file is in:  #python3 -m src <filename of a desired python module file without the .py>

The inputs folder contains the input files where the user can update the inputs the program uses.

the analysis_mode list is used to select the types of models used by the program.

[TODO: DOCUMENT MODELS IN EACH FILE]

source venv/bin/activate


## Context:
This is part of my personal rocket engine design project. It started as a continuation of Ben Klammer's work [1]
He does some really cool things like optimization that this program does not do. 
I would recommend checking his thesis and program out. His work taught me a good amount of what I know about propulsion

I decided to make this program over already built software packages because my small scale engine's design loads are very dependent on the self-pressurized
feed system, and I want to build design tools that will help me understand performance when I build my engine.

Through working with different system models, the program has evolved into a relatively modular framework where the user can easily add new models.

I have definitely learned a lot through making this program. Even when the program is broken, its cool to see how the errors propagate through the program. I have read about vehicle failures where the failure of a component causes several effects that propagate through the system, and I definitely see that when I have extraneous inputs.

I got some of the liquid stuff to work and accidentally broke the hybrid models. I likely won't fix the hybrid stuff for a while.

#Use at your own risk lol, I am still developing and debugging this program. When things are more stable I will do better at documentation


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
| [19]  | Modeling Feed System Flow Physics for Self-Pressurizing Propellants, Stanford | Span Wagner Coeffs | https://www.researchgate.net/publication/268482381_Modeling_Feed_System_Flow_Physics_for_Self-Pressurizing_Propellants |
| [20]  | AIAA Liquid Rocket Thrust Chambers | | |
| [21]  | Huzel + Huang | | |
| [22]  | Ben Klammer 446 CC model | | |
| [23]  | Implementing CEA calculations using Cantera - Kyle Niemeyer | Starting to outgrow rocketcea, looking into using cantera for cc model because it is more compatible with differential eqns | https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html |
| [24]  | NASA SP-1311-2 | fixing paraffin definition in rocketcea. see PARAFFIN_DEFFINITION.md | https://shepherd.caltech.edu/EDL/PublicResources/sdt/refs/NASA-RP-1311-2.pdf |



Ryan Wright, From Calgary, Alberta, Canada   

<!--pyenv shell king-hybrid-env -->