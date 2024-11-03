### Rocket Engine Simulation 🚀

## RUNNING THE PROGRAM:

In the directory this file is in:  #python3 -m src <filename without the .py>

The inputs folder contains files where the user can update the inputs the program uses.

the analysis_mode vector selects the types of models used by the program.
[TODO: DOCUMENT MODELS IN EACH FILE]



## Context:
A continuation of Ben Klammer's work. All the hybrid stuff in this program is ripped from of it:
https://github.com/UVicRocketry/HybridModeling
He does some really cool things like optimization that this program does not do. 
I would recommend checking his thesis and program out. His work taught me a good amount of what I know about propulsion

I made this progam to help with rocket engine design. I decided to make this over alternatives like
already built software packages alternatives because I enjoy programming and I value having the control
of creating my own models and the understanding it brings me. I think it is important to have a good amount
of understanding and control because of how finicky nitrous oxide is.

Maybe this isn't the most efficient choice though. I have definitely learned a lot through making it though. Even when the program is broken, its cool to see how the errors propagate through the program. I have read about vehicle failures where the failure of a component causes several effects that propagate through the system, and I definitely see that when I have extraneous inputs.

This program was designed to be a framework for a user to combine and use different models to simulate a rocket engine.

I got some of the liquid stuff to work and by doing that broke the hybrid models. I am no longer interested in hybrids so I likely won't fix the hybrid stuff.

###Use at your own risk lol, I am still developing this program 

## Sources and Citations:
| Number | Source                                                                                                                                                         |
|--------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1      | Ben Klammer, Hybrid Modelling: https://github.com/bklammer/HybridModeling                                                                                      |
| 2      | Emerson Vargas Niño, Mohammad Reza H. Razavi, Design of Two-Phase Injectors Using Analytical and Numerical Methods with Application to Hybrid Rockets https://emersonvn.com/project/two_phase_injector/# |
| 3      | Benjamin S. Waxman, Jonah E. Zimmerman, Brian J. Cantwell, Mass Flow Rate and Isolation Characteristics of Injectors for Use with Self-Pressurizing Oxidizers in Hybrid Rockets https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf |
| 4      | Zilliac, Gregory & Karabeyoglu, Arif. (2005). Modeling of Propellant Tank Pressurization.                                                                      |



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
