### Rocket Engine Simulation 🚀

This program was designed to be a framework for a user to combine and use different models to simulate a rocket engine.

## RUNNING THE PROGRAM:

In the directory this file is in:  #python3 -m src <filename without the .py>

The inputs folder contains the input files where the user can update the inputs the program uses.

the analysis_mode vector selects the types of models used by the program.
[TODO: DOCUMENT MODELS IN EACH FILE]



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
| Number   | Source                                                                                                                                                         |
|----------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [1]      | Ben Klammer, Hybrid Modelling: https://github.com/bklammer/HybridModeling                                                                                      |
| [2]      | Emerson Vargas Niño, Mohammad Reza H. Razavi, Design of Two-Phase Injectors Using Analytical and Numerical Methods with Application to Hybrid Rockets https://emersonvn.com/project/two_phase_injector/# |
| [3]      | Benjamin S. Waxman, Jonah E. Zimmerman, Brian J. Cantwell, Mass Flow Rate and Isolation Characteristics of Injectors for Use with Self-Pressurizing Oxidizers in Hybrid Rockets https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf |
| [4]      | Zilliac, Gregory & Karabeyoglu, Arif. (2005). Modeling of Propellant Tank Pressurization.                                                                      |



Ryan "king hybrid" Wright, From Calgary, Alberta, Canada --> University of Victoria, BC, Canada, 2023 --> forever (mech lab reference)

```plaintext
                                             &@&.                       
                                          @@      /@                    
                               %@@@@@@,  @&    @%   %(                  
                           (@%         @@@        @                     
              ,&@@@@@@@@@@@.         @@&         @#                     
          *@@@@@@&      @/         @@,       ,&,  /@@@.                 
         @@@@@%        @    &@@@@@@.                 @@%                
        #@@@@@        @..@*    @@                     @@                
        *@@@@@        @,    (@/                      &@,                
         @@@@@@          @@.         *@@@@@,        #@#                 
          @@@@@@    (@#           #@@      @       @@.                  
            @@@@@@  .&@@@@@@    @@ @      @/     /@&                    
             #@@@@@@.    #@   &@  @      @     @@/  #@,                 
               .@@@@@@@. @@  @@@  @    @.   @@%     @@@%                
               @  @@@@@@@@@ % @  ,   @%@@@*         #@@@                
              @      %@@@@@@@@@.                    @@@@/                       
             @           @@@@@@@@@@@@,           (@@@@@@                
             @          *@.  *@@@@@@@@@@@@@@@@@@@@@@@@@                 
            @/      .@@            ,@@@@@@@@@@@@@@@@                    
           @    @@,                                                     
          @@%                                                                              
