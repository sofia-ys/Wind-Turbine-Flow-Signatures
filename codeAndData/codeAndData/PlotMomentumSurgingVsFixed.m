clear all
warning('off', 'all');
folder = "C:\Users\janwa\Documents\TU Delft\Bsc2\Test simulation\Code\codeAndData\NREL_FXXXXX_5D_000_00025_copy\exportedCSV\phase_average\";
planes = [2, 5, 8];
momentum_case = 'Fixed';

momentum_fixed = zeros(3, 24);
instances = 24;
for i = 1:length(planes)
   file = (folder +'PhaAvg_V_comp_' + momentum_case+ '_xPlane'+num2str(planes(i))+"D.csv");
   data = readtable(file);
   if i == 1
       col1 = data{:, end-1};  % Second to last column
       col2 = data{:, end};    % Last column

       distance = sqrt(col1.^2 + col2.^2);
       rows_to_keep = find(distance <= 126);
   end
   for ii = 1:instances
       if ii == 1 & i == 1
            data = data(rows_to_keep, :);
       end
       momentum_fixed(i,ii) = sum(data{:,(3*ii-2)});
   end
end
momentum_surging = zeros(3, 24);
momentum_case = 'Surging';
for i = 1:length(planes)
   file = (folder +'PhaAvg_V_comp_' + momentum_case+ '_xPlane'+num2str(planes(i))+"D.csv");
   data = readtable(file);
   if i == 1
       col1 = data{:, end-1};  % Second to last column
       col2 = data{:, end};    % Last column

       distance = sqrt(col1.^2 + col2.^2);
       rows_to_keep = find(distance <= 126);
   end
   for ii = 1:instances
       if ii == 1 & i == 1
            data = data(rows_to_keep, :);
       end
       momentum_surging(i,ii) = sum(data{:,(3*ii-2)});
   end
end

% Example arrays - mock data 
fixedData = momentum_fixed;   % Example random data for fixed case (3 planes, 24 phases)
surgingData = momentum_surging; % Example random data for surging case (3 planes, 24 phases)

% Phase angles from 0 to 2*pi (24 phases)
phase = linspace(0, 2*pi, 24);  % 24 values from 0 to 2*pi

% Initialize the figure
figure;

% Loop over the 3 planes (2D, 5D, and 8D)
for i = 1:3
    subplot(3,1,i);
    
    plot(phase, fixedData(i,:), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    hold on;
    plot(phase, surgingData(i,:), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    
    title(['Momentum vs Phase for Plane ' num2str(i*3-1) 'D']);
    xlabel('Phase (rad)');
    ylabel('Momentum');
    legend('Fixed Case', 'Surging Case');

    % Set x-axis limits and ticks
    xlim([0, 2*pi]);
    xticks([0 pi/2 pi 3*pi/2 2*pi]);
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
    
    hold off;
end
