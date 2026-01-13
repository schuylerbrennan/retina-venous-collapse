% To run the compartmental model code, call
% FullCompartmentSystem. The main components to be changed in this file are M0_vary and
% Pa_vary, which can be set to true or false. M0v and Pav are two vectors
% that can be changed to contain the values of M0 (oxygen saturation) and Pa
% (incoming arterial pressure)you would like to include.

% The vector IOPv should contain a list of the IOP values you would
% like to run. The IOP values should be in dyn/cm^2.

% The variable colors can be changed to contain a string of color codes
% indicating the color of the plot curves. When plotting, the order of the color
% codes will follow the order of the IOP values in IOPv.

% To turn venous collapsibility on and off, you will need to go to
% obtainnewflowquantselev. To turn collapsibility on, uncomment the if
% statements in lines 147-170 and comment out lines 174-180. To turn
% collapsibility off, comment out lines 147-170 and uncomment lines
% 174-180. 

% In postprocessorPa and postprocessorM0, there are for loops that are
% commented out that will prevent the code from plotting trials with
% negative transmural pressure in the capillaries past a certain threshold.
% To use this feature, uncomment the for loop, change the input to the
% function to read Pa1 = varargin{3} or M01 = varargin{3}, and change to n
% = numel(Pa1) or n = numel(M01).


