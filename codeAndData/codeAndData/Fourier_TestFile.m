% CLEAR PREVIOUS WORKSPACE
clear; close all; clc;

% DEFINE PARAMETERS
folder = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_y00D";
files = dir(fullfile(folder, '*.csv')); % Get all CSV files
numFiles = length(files);  % Adjusted for skipped files

% Ensure we have valid files
if numFiles < 1
    error("No CSV files found in the folder!");
end

% Read the first file to determine grid size
firstFilePath = fullfile(folder, files(1).name);
firstData = readmatrix(firstFilePath);

if isempty(firstData) || size(firstData, 2) < 5
    error("First file is empty or has insufficient columns!");
end

numGridPoints = size(firstData, 1);  % Rows represent grid points

% Preallocate velocity structure array (each grid point has its velocity and metadata)
V_struct = struct('velocity', nan(numFiles, 1), 'metadata', cell(numFiles, 1));  % Preallocate for numFiles time steps

% Start parallel pool if not already open
pool = gcp('nocreate'); 
if isempty(pool)
    parpool; % Automatically selects best number of workers
end

% Start timing
totalStartTime = tic; 
fprintf("Processing %d files with %d grid points...\n", numFiles, numGridPoints);

% LOOP THROUGH FILES (TIME STEPS) using parfor
parfor t = 1:numFiles
    filePath = fullfile(folder, files(t).name);
    
    % Use readmatrix for faster reading
    data = readmatrix(filePath);

    % Ensure data is valid
    if isempty(data) || size(data, 2) < 5
        continue;  % Skip invalid files
    end

    % Extract velocity components
    u_velocity = data(:, 3);
    v_velocity = data(:, 4);
    w_velocity = data(:, 5);

    % Compute absolute velocity magnitude
    velocity_magnitude = sqrt(u_velocity.^2 + v_velocity.^2 + w_velocity.^2);

    % Store velocity data in structure array
    V_struct(t).velocity = velocity_magnitude;  % Store velocity for current time step

    % Optionally, add metadata related to the file (e.g., filename, time step)
    V_struct(t).metadata = struct('file_name', files(t).name, 'time_step', t * 0.013778915 * 15);
    
    % Progress update every 100 files
    if mod(t, 100) == 0
        fprintf("Processed %d/%d time steps...\n", t, numFiles);
    end
end

% Display total elapsed time
totalElapsedTime = toc(totalStartTime);
fprintf("Processing complete! Total time: %.2f seconds.\n", totalElapsedTime);

% Save data in .mat format
save('V_velocity_struct.mat', 'V_struct', '-v7.3'); % v7.3 enables compression for large data
fprintf("Data saved to 'V_velocity_struct.mat'.\n");
