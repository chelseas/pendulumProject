README
Chelsea Sidrane
4730 Final Project

To use the simulation:
1. Open the script FPwrite.m in MATLAB
2. Type ‘close all; clear all; clc;’at the command window
3. Run FPwrite by clicking the green run button. FPwrite takes no arguments.

IMPORTANT NOTE
Before running FPwrite and between EACH RUN of FPwrite, you must type ‘close all; clear all; clc;’  

Detailed Description
FPwrite will run a demonstration of the suite of code. The triple pendulum is first solved three ways and all animations are simulated simultaneously. Then a limiting case is demonstrated where the mass of the end link is very large and thus the motion approximates simple pendulum motion. The four bar linkages are then demonstrated. The motion is first simulated for very simple initial conditions where the four bar linkage swings in a predictable fashion. The motion of the four bar linkage is then simulated using randomized lengths, randomized masses, and randomized initial angles. Lastly, an n-pendulum is simulated with 10 links. 

Notes for the Curious
Feel free to change the initial conditions specified for the simulations prescribed in FPwrite. Also feel free to remove the ‘pause()’ commands between plots popping up. These were added for demonstration purposes because 19 figures popping up all once is annoying. Also also feel free to change the number of links in the n-pendulum simulation in FPwrite.