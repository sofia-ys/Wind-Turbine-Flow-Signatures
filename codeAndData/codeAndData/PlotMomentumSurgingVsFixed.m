% Example arrays - mock data 
fixedData = rand(3, 24);   % Example random data for fixed case (3 planes, 24 phases)
surgingData = rand(3, 24); % Example random data for surging case (3 planes, 24 phases)

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