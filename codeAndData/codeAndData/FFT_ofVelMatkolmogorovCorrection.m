clear all
load('V_velocity_matrix.mat');  

velocityData = V_matrix.velocity;  % Access the velocity field

[gridPoints, timeSteps] = size(velocityData);

% Initialize the array to store the dominant frequency for each grid point
dominantFrequencies = zeros(gridPoints, 1);

dt = 0.013778915 * 15;  % based on simulation data

% Calculate the Fourier transform for each grid point (row)
for i = 1:gridPoints
    velocity = velocityData(i, :); % selecting the index for the row and selecting all columns
    
    % Apply the Fourier Transform (fft) to the velocity data
    freqDomain = fft(velocity);
    
    % Calculate the corresponding frequencies
    frequencies = (0:timeSteps-1) / (timeSteps * dt);
    
    % Get the magnitude of the frequency spectrum
    magnitude = abs(freqDomain);
    
    % Apply Kolmogorov correction: compensate by f^(5/3)
    % Avoid f=0 (DC component), so skip the first element
    kolmogorovCorrection = frequencies .^ (5/3);
    kolmogorovCorrection(1) = 0; % prevent Inf or NaN for DC component

    adjustedMagnitude = magnitude .* kolmogorovCorrection;
    
    % Find the index of the maximum adjusted magnitude (dominant frequency)
    [~, maxIndex] = max(adjustedMagnitude(2:end));  % Skip DC component (first element)
    
    % Store the dominant frequency
    dominantFrequencies(i) = frequencies(maxIndex + 1);  % +1 to adjust index
end

% Save the dominant frequencies as a .mat file
save('dominantFrequencies.mat', 'dominantFrequencies');
