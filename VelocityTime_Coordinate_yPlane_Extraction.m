% CLEAR PREVIOUS WORKSPACE
clear; close all; clc;

% DEFINE PARAMETERS
folder = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_y00D";
files = dir(fullfile(folder, '*.csv')); % Get all CSV files
numFiles = length(files);  

% Ensure we have valid files
if numFiles < 1
    error("No CSV files found in the folder!");
end

% User-defined coordinate columns
x_column = 9;  
z_column = 11;  

% Read first file to determine grid size
firstFilePath = fullfile(folder, files(1).name);
firstData = readmatrix(firstFilePath, 'Delimiter', ','); % Specify delimiter for speed

if isempty(firstData) || size(firstData, 2) < 5
    error("First file is empty or has insufficient columns!");
end

numGridPoints = size(firstData, 1);  

% Preallocate a numeric matrix for velocity storage
velocity_data = nan(numGridPoints, numFiles, 'single'); 

% Read x and z coordinates once from the first file
x_coords = firstData(:, x_column);
z_coords = firstData(:, z_column);

% Start parallel pool if not already open
numWorkers = 15;  
pool = gcp('nocreate'); 
if isempty(pool)
    parpool(numWorkers); 
end

% Start timing
totalStartTime = tic; 
fprintf("Processing %d files with %d grid points...\n", numFiles, numGridPoints);

% LOOP THROUGH FILES (TIME STEPS) using parfor
parfor t = 1:numFiles
    filePath = fullfile(folder, files(t).name);
    
    % Use readmatrix for faster reading
    data = readmatrix(filePath, 'Delimiter', ',');

    % Ensure data is valid
    if isempty(data) || size(data, 2) < 5
        continue;  
    end

    % Compute absolute velocity magnitude
    velocity_data(:, t) = sqrt(data(:, 3).^2 + data(:, 4).^2 + data(:, 5).^2);

    % Progress update every 100 files
    if mod(t, 100) == 0
        fprintf("Processed %d/%d time steps...\n", t, numFiles);
    end
end

% Create a structure with metadata
V_matrix = struct();
V_matrix.velocity = velocity_data;  % 2D matrix: [grid points x time steps]
V_matrix.x = x_coords;              % x-coordinates (1D array)
V_matrix.z = z_coords;              % z-coordinates (1D array)
V_matrix.time_steps = numFiles;      % Number of time steps
V_matrix.grid_points = numGridPoints; % Number of grid points

% Display total elapsed time
totalElapsedTime = toc(totalStartTime);
fprintf("Processing complete! Total time: %.2f seconds.\n", totalElapsedTime);
fprintf("SAVING FILE. PLEASE WAIT FOR THIS PROCESS TO COMPLETE!!!\n");

% Save using efficient format
save('V_velocity_matrix.mat', 'V_matrix', '-v7.3');  
fprintf("Data saved to 'V_velocity_matrix.mat'.\n");
