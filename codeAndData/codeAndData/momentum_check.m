clear all
close all
warning('off', 'all');
folder_base = "C:\Users\janwa\Documents\TU Delft\Bsc2\Test simulation\Code\codeAndData\NREL_FXXXXX_5D_000_00025_copy\exportedCSV\";
planes = [2, 5, 8];
momentum_case = 'Fixed';
radius_length = 63;
instances = 2400;
momentum_surging = zeros(3, instances);
%energy_fixed = zeros(3, 24);

for i = 1:length(planes)
   folder = (folder_base +'exportedCSV_x_'+num2str(planes(i))+"D\");
   for ii = 1:instances
       file = (folder+ 'xPlane_'+num2str(planes(i))+'D_'+ sprintf('%04d', ii)+'.csv');
       data = readtable(file);
       col1 = data{:, end-1};  % Second to last column
       col2 = data{:, end};    % Last column

       distance = sqrt(col1.^2 + col2.^2);
       rows_to_keep = find(distance <= radius_length);
       data = data(rows_to_keep, :);
       momentum_surging(i,ii) = mean(data{:,(3)});
       if mod(ii, 100) == 0
           disp(ii)
       end
   end
end
save(['mean_velocity_', momentum_case,'.mat'], 'momentum_surging'); 
time = zeros(1,leninstances);
for i = 0:(len(instances)-1)
    time(1,i) = i*0.013778915 * 15;
end
figure;
plot(momentum_surging(1,:),  '-r', 'LineWidth', 1.5); hold on;
plot(momentum_surging(2,:), '-g', 'LineWidth', 1.5);
plot(momentum_surging(3,:), '-b', 'LineWidth', 1.5);
hold off;

xlabel('Time');
ylabel('Mean velocity');
title('Mean Velocity Over Time');
legend('2D', '5D', '8D');
grid on;
exportgraphics(FigMaster,"./NREL_FXXXXX_5D_000_00025_copy/exportedPNG/mean_velocity_surging.png")
