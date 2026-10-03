function filtered_data = lowpass_filter(data, fs, fc)
    % lowpass_filter applies a Butterworth low-pass filter to the input data.
    %
    % Inputs:
    %   data - An N-dimensional array with time along the first dimension.
    %   fs   - Sampling rate (Hz)
    %   fc   - Cutoff frequency (Hz)
    %
    % Output:
    %   filtered_data - The low-pass filtered data.
    
    % Filter design parameters
    order = 4;                 % Filter order; adjust as needed
    Wn = fc/(fs/2);            % Normalized cutoff frequency
    
    % Create Butterworth filter coefficients for a low-pass filter
    [b, a] = butter(order, Wn, 'low');
    
    % Apply zero-phase filtering along the first dimension
    filtered_data = filtfilt(b, a, data);

	if false
		%% Demo for low pass filtering using a Butterworth filter
		% Sampling parameters
		fs = 1000;         % Sampling rate in Hz
		T  = 1;            % Duration in seconds
		t  = 0:1/fs:T-1/fs; % Time vector

		% Create sample data:
		% A 5 Hz low-frequency signal plus a 100 Hz high-frequency noise component.
		signal  = sin(2*pi*5*t);
		noise   = 0.5*sin(2*pi*100*t);
		data    = signal + noise;
		data    = data(:); % Ensure data is a column vector (time along first dimension)

		% Define cutoff frequency for low pass filtering
		fc = 10;  % Cutoff frequency in Hz (choose to remove 100 Hz noise)

		% Apply low pass filtering using the custom function
		filtered_data = lowpass_filter(data, fs, fc);

		% Plot the results
		figure;
		plot(t, data, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Original Signal');
		hold on;
		plot(t, filtered_data, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Filtered Signal');
		xlabel('Time (s)');
		ylabel('Amplitude');
		title('Low Pass Filter Demonstration');
		legend;
		grid on;
	end
end
