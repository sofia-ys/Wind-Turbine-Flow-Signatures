% CLEAR UP PREVIOUS WORKSPACE

clear all 
close all


% DEFINE SOME GENERAL PARAMETERS
timeStepSize = 0.013778915 * 15;

tic

saveFigure = 1;

closeAllIndex = 0;

D = 126;
V0 = 11.4;

% switch for if always plot and save
alwayPlotting = 1;


% DEFINE PLOTTING COORDINATE SPACE

% Define the stepsize
stepsize = 1/80;

% WARNING: Coordinates are defined in x/D !!!!

% Defining the boundaries of the plotting domain
xMin = -2 + stepsize/2;
xMax = 8 + stepsize/2;

zMin = -1 + stepsize/2;
zMax = 1 + stepsize - stepsize/2;

[x_grid,z_grid] = meshgrid( ...
    xMin:stepsize:xMax, ...
    zMin:stepsize:zMax);

% Establishing the mesh matrices for the coordinates in x and z
coordinate_yPlane.x_grid = x_grid;
coordinate_yPlane.z_grid = z_grid;





% LOCATE THE DATA FILES 

folder = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_y00D";

%/yPlane_" + indexOfCSVfileString + ".csv"

files = dir(fullfile(folder, '*.csv')); % Stores all the names of each csv file in the folder


% DEFINE DATA STORAGE ARRAYS
V_velocity = zeros(length(files), 1); % Establish empty array to store absolute velocity magnitude 


for i = 1:length(files) % For each csv file

    % Open file(i) for reading
    file = fopen(fullfile(folder, files(i).name), 'r');

    % Skip first 9 rows
    for j = 1:1
        fgetl(file);
    end

    % Read ONLY row 10 (ignore everything else)
    rowData = fgetl(file); % Read row 10 as a single string

    % Close the file immediately after reading row 10
    fclose(file);

    % Convert rowData string to numbers and extract column 5
    numericData = str2double(strsplit(rowData, ',')); % Convert CSV string to numeric array

    % Extract the value of each velocity component
    u_velocity = numericData(3); % Extract the velocity of u-component from column 3
    v_velocity = numericData(4);
    w_velocity = numericData(5);
    
    % Compute the absolute velocity magnitude and store it in the array
    V_velocity(i) = sqrt(u_velocity^2 + v_velocity^2 + w_velocity^2);  
end



%velocities = cell2mat(velocities);

T = timeStepSize; % Sampling period
Fs = 1/T; % Sampling frequency (Hz), change based on your data
N = length(V_velocity); % Number of data points

t = (0:N-1) * T; % Time vector

% Compute FFT
U = fft(V_velocity); % Perform FFT
f = Fs * (0:(N/2)) / N; % Frequency vector (one-sided spectrum)
P2 = abs(U/N); % Two-sided spectrum
P1 = P2(1:N/2+1); % Single-sided spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Adjust magnitude

[peakSpectrum, peakFrequencies] = findpeaks(f); % Find local maxima

% Plot time-domain signal
figure;
subplot(2,1,1);
plot(t, V_velocity);
title('Time-Domain Signal');
xlabel('Time (s)');
ylabel('Absolute Velocity');
grid on;

subplot(2,1,2);
semilogy(f, P1); % Use logarithmic scale for y-axis
title('Frequency-Domain (FFT) - Log Scale');
xlabel('Frequency (Hz)');
ylabel('|U(f)| (Log Scale)');
grid on;