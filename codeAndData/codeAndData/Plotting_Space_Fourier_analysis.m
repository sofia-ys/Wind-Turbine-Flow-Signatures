% THIS CODE DEFINES THE PLOTTING COORDINATE SPACE FOR THE PROMINENT
% FREQUENCIES 

% Define the step size in the plotting coordinate space
stepsize = 1/80;

% WARNING COORDINATES ARE DEFINED IN TERMS OF X/D !!!

% Defining the boundaries of the plotting domain
x_min = -2 + stepsize/2;
x_max = 8 + stepsize/2;

z_min = -1 + stepsize/2;
z_max = 1 + stepsize - stepsize/2;

% Establishing the mesh matrices for the coordinates in x and z
[ x_grid , z_grid ] = meshgrid( ...
    x_min : stepsize : x_max, ...
    z_min ...
    : stepsize : z_max);


coordinate_yPlane.x_grid = x_grid;
coordinate_yPlane.z_grid = z_grid;
