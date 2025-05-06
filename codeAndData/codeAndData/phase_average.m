prompt = 'Filename of plane you want to phase-average: ';
input_filename = input(prompt, 's');
load(['xPlane Velocity component field vs time\' input_filename '.mat']);

fprintf("Running...\n")

% Amount of coordinate points
num_points = size(V_matrix.velocity, 1);
fprintf("Read amount of coordinate points: %d\n", num_points)

% Saving flags
save_matfile = 1;
save_csvfile = 1;

% Phase averaging: overwrite first 24 time frames
for coordinate_point = 1:num_points
    if mod(coordinate_point, 10000) == 0
        fprintf("Processing coordinate point %d / %d\n", coordinate_point, num_points)
    end
    for component = 1:3
        for phase_time = 1:24
            comp_sum = 0;
            for period = 0:99 
                comp_sum = comp_sum + V_matrix.velocity(coordinate_point, component, phase_time + period * 24);
            end
            comp_avg = comp_sum / 100;
            V_matrix.velocity(coordinate_point, component, phase_time) = comp_avg;
        end
    end
end

% Truncate the velocity data to just the 24 averaged frames
V_matrix.velocity = V_matrix.velocity(:, :, 1:24);

% Save the reduced MAT file
if save_matfile == 1
    output_filename = sprintf('PhaAvg_%s.mat', input_filename);
    save(output_filename, 'V_matrix', '-v7.3');
    fprintf("Data saved to '%s'.\n", output_filename);
end

% Save to CSV (velocity + coordinates, no headers)
if save_csvfile == 1
    velocity_2D = reshape(V_matrix.velocity, num_points, []);
    
    % Use correct coordinates depending on plane orientation
    y = V_matrix.y;
    z = V_matrix.z;
    
    full_data = [velocity_2D, y, z];x
    
    csv_filename = sprintf('PhaAvg_%s.csv', input_filename);
    writematrix(full_data, csv_filename);
    
    fprintf("Data saved to CSV (no headers): %s\n", csv_filename);
end
