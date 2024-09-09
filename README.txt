### RUNNING THE PROGRAM:

python3 -m src <filename without the .py> in the directory this file is in. 

The inputs folder contains files where the user can update the inputs the program uses.

the analysis_mode vector selects the types of models used by the program.
[TODO: DOCUMENT MODELS IN EACH FILE]



### Context:
A continuation of Ben Klammer's work. All the hybrid stuff in this program is based off of it:
https://github.com/UVicRocketry/HybridModeling
He does some really cool things like optimization that this program does not do. 
I would recommend checking his thesis and program out.

I made this progam to help with rocket engine design. I decided to make this over alternatives like
already built software packages alternatives because I enjoy programming and I value having the control
of creating my own models and the understanding it brings me. 

Maybe this isn't the most efficient choice tho.

This program was designed to be a framework for a user to combine and use different models to simulate a rocket engine.
[In Progress] - Users can add models to the models folder and update the map in the thrust curve script.
[untested/unfinished] - liquid rocket engine models
[sensitivity analysis is not supported for non hybrid yet]
Program has tested and validated hybrid component models and in progress liquid models




%                                             &@&.                       
%                                          @@      /@                    
%                               %@@@@@@,  @&    @%   %(                  
%                           (@%         @@@        @                     
%              ,&@@@@@@@@@@@.         @@&         @#                     
%          *@@@@@@&      @/         @@,       ,&,  /@@@.                 
%         @@@@@%        @    &@@@@@@.                 @@%                
%        #@@@@@        @..@*    @@                     @@                
%        *@@@@@        @,    (@/                      &@,                
%         @@@@@@          @@.         *@@@@@,        #@#                 
%          @@@@@@    (@#           #@@      @       @@.                  
%            @@@@@@  .&@@@@@@    @@ @      @/     /@&                    
%             #@@@@@@.    #@   &@  @      @     @@/  #@,                 
%               .@@@@@@@. @@  @@@  @    @.   @@%     @@@%                
%               @  @@@@@@@@@ % @  ,   @%@@@*         #@@@                
%             /#      %@@@@@@@@@.                    @@@@/                       
%            /%           @@@@@@@@@@@@,           (@@@@@@                
%             @          *@.  *@@@@@@@@@@@@@@@@@@@@@@@@@                 
%            @/      .@@            ,&@@@@@@@@@@@@@@@                    
%           @    @@,                                                     
%          @@%                                          