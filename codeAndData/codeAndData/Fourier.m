clear all 
close all

timeStepSize = 0.013778915 * 15;

tic

saveFigure = 1;

closeAllIndex = 0;

% Define the stepsize in the Coordinate Space
stepsize = 1/80;


% WARNING COORDINATES ARE DEFINED IN TERMS OF X/D !!!

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

%%  Which cases to compare, basic settings

D = 126;
V0 = 11.4;


figureFileNameMaster = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/";

fileNumber = 1;

figNumber = fileNumber;

%% interp data to uniform grid
how_many = 1600;
% indexOfCSVfileArray = 0:1:2400;
indexOfCSVfileArray = 0:1:how_many;

% switch for if always plot and save
alwayPlotting = 1;

%% Read the data
u_1 = cell(how_many+1, 1);
for indexCSVHolder = indexOfCSVfileArray

    indexOfCSVfile = indexCSVHolder;
    indexOfCSVfileString = sprintf( '%04d', indexOfCSVfile );

    csvTail = {};
    csvTail{1} = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_y00D/yPlane_" + indexOfCSVfileString + ".csv";
    numberOfPlanes = length(csvTail);
    fileName_All = cell(fileNumber,numberOfPlanes);

    for pp = 1:fileNumber
        for gg = 1:numberOfPlanes
            fileName_All{pp, gg} = csvTail{gg};
        end
    end
    
     %% Read in all the needed data

    dataIn = cell(fileNumber,numberOfPlanes);


    for tt = 1:fileNumber

        for gg = 1:numberOfPlanes
            dataIn{tt, gg} = readtable(fileName_All{tt, gg}, 'VariableNamingRule', 'preserve');
        end

    end


    openFaomU = cell(fileNumber,numberOfPlanes);
    openFaomV = cell(fileNumber,numberOfPlanes);
    openFaomW = cell(fileNumber,numberOfPlanes);


    openFaomX = cell(fileNumber,numberOfPlanes);
    openFaomY = cell(fileNumber,numberOfPlanes);
    openFaomZ = cell(fileNumber,numberOfPlanes);
    
    for pp = 1:fileNumber

        for gg = 1:numberOfPlanes

            openFaomU{pp, gg} = dataIn{pp, gg}.("U:0");
            openFaomV{pp, gg} = dataIn{pp, gg}.("U:1");
            openFaomW{pp, gg} = dataIn{pp, gg}.("U:2");

            openFaomX{pp, gg} = dataIn{pp, gg}.("Points:0");
            openFaomY{pp, gg} = dataIn{pp, gg}.("Points:1");
            openFaomZ{pp, gg} = dataIn{pp, gg}.("Points:2");
            
            
            u_1{indexCSVHolder+1} = openFaomU{pp, gg}(1);
         

        end

    
    end

    
    clear dataIn;
end
u_1 = cell2mat(u_1);

T = timeStepSize; % Sampling period
Fs = 1/T; % Sampling frequency (Hz), change based on your data
N = length(u_1); % Number of data points

t = (0:N-1) * T; % Time vector

% Compute FFT
U = fft(u_1); % Perform FFT
f = Fs * (0:(N/2)) / N; % Frequency vector (one-sided spectrum)
P2 = abs(U/N); % Two-sided spectrum
P1 = P2(1:N/2+1); % Single-sided spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Adjust magnitude

% Plot time-domain signal
figure;
subplot(2,1,1);
plot(t, u_1);
title('Time-Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
semilogy(f, P1); % Use logarithmic scale for y-axis
title('Frequency-Domain (FFT) - Log Scale');
xlabel('Frequency (Hz)');
ylabel('|U(f)| (Log Scale)');
grid on;




