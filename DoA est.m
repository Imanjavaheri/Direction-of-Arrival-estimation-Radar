%% Part 3-MATLAB Implementation
clc; clear; close all;

% Define parameters
N = 58; % Number of antennas
lambda = 3.896e-3; % Wavelength 

dx = lambda/2; % Antenna spacing 
% dx = lambda;
antenna_positions = (0:N-1) * dx;  

k = 2 * pi / lambda; % Wavenumber

%Single target
% theta_true = 10;

%multiple target
theta_true = [15, 17]; % Angles of two targets in degrees 

theta_rad = deg2rad(theta_true); % Convert to radians
K = length(theta_true); % Number of targets

% Reshape vectors to align dimensions for broadcasting
antenna_positions = antenna_positions(:); % N x 1 column vector
theta_rad = theta_rad(:)'; % 1 x K row vector

% Compute phase shifts matrix
% phase_shifts N x K matrix
phase_shifts = -k * antenna_positions * sin(theta_rad); % N x K matrix

% Compute received signals from all targets
% Each column corresponds to a target
received_signals = exp(1j * phase_shifts); % N x K matrix

% Sum the signals across targets for each antenna
received_signal = sum(received_signals, 2); % N x 1 vector

% Add noise
sigma_noise = 0.1;
noise = sqrt(sigma_noise/2) * (randn(N, 1) + 1j * randn(N,1));
received_signal = received_signal + noise;


% Define scanning angles
theta_scan = -90:0.1:90; % degrees
theta_scan_rad = deg2rad(theta_scan);  

% Initialize array to store beamforming output
beamforming_output = zeros(length(theta_scan), 1);

% Loop over scanning angles
for idx = 1:length(theta_scan_rad)
    % Compute steering vector for the scanning angle
    steering_vector = exp(-1j * k * antenna_positions * sin(theta_scan_rad(idx)));
    
    % Compute beamformer output (power)
    beamforming_output(idx) = abs(steering_vector' * received_signal)^2;
end

% Normalize the beamforming output
beamforming_output = beamforming_output / max(beamforming_output);

% Find peaks in the beamforming output
[pks, locs] = findpeaks(beamforming_output, 'MinPeakHeight', 0.5, 'SortStr', 'descend', 'NPeaks', K);

% Estimated angles
theta_est = theta_scan(locs);

% Sort the estimated angles
theta_est = sort(theta_est);

% Plot the Beamforming Spectrum

figure;
plot(theta_scan, beamforming_output, 'LineWidth', 1.5);
xlabel('Angle (degrees)');
ylabel('Normalized Power');
title('Beamforming Output');
grid on;
hold on;
plot(theta_est, pks, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Beamforming Output', 'Estimated Angles');

% Output Results
fprintf('True Angles: %.2f degrees\n', theta_true);
fprintf('Estimated Angles: %.2f degrees\n', theta_est);


%% part 4-2D Position Estimation
clc; clear; close all;

% Parameters
N = 58;
lambda = 3.896e-3;
dx = lambda/2;  
antenna_positions = (0:N-1) * dx; 

% Target Parameters
theta_true = 26;  
theta_rad = deg2rad(theta_true); 
r = 83; % Target range in meters

% Received Signal Parameters
k = 2 * pi / lambda; % Wavenumber

% Time Delay Calculation
tau = 2 * r / 3e8; % Round-trip time delay (seconds)

% Sampling Parameters
B = 1e9; % Bandwidth
Ts = 1 / (10 * B); % Sampling interval (oversampling by factor of 10)
Fs = 1 / Ts; % Sampling frequency

% Time Vector 
t_max = tau + 1e-6;  
t = 0 : Ts : t_max;  

% Generate Transmitted Signal
g_tx = sinc(B * t) .* exp(1j * 2 * pi * 77e9 * t);  

% Simulate Received Signal at Each Antenna

% Initialize received signal matrix (N antennas x length of time vector)
received_signal = zeros(N, length(t));

% Correct Phase Shifts
phase_shifts = -k * antenna_positions * sin(theta_rad); % 1 x N vector

% Generate Received Signals
for n = 1:N
    % Time-delayed signal
    t_delayed = t - tau;
    % Initialize signal with zeros
    signal = zeros(1, length(t));
    % Identify indices where t_delayed >= 0
    idx_positive = find(t_delayed >= 0);
    % Assign the delayed signal 
    signal(idx_positive) = sinc(B * t_delayed(idx_positive)) .* ...
                            exp(1j * (2 * pi * 77e9 * t_delayed(idx_positive) + phase_shifts(n)));
    % Received signal at antenna n
    received_signal(n, :) = signal;
end

% Add Noise
SNR_dB = 30; 
noise_power = 10^(-SNR_dB/10);
noise = sqrt(noise_power/2) * (randn(N, length(t)) + 1j * randn(N, length(t)));
received_signal_noisy = received_signal + noise;

% Demodulate Received Signals to Baseband

received_signal_demod = received_signal_noisy .* repmat(exp(-1j * 2 * pi * 77e9 * t), N, 1);

% Perform Matched Filtering for Range Estimation
ref_signal = sinc(B * t);

% Initialize range profile matrix
range_profile = zeros(N, length(t));

% Perform matched filtering using cross-correlation
for n = 1:N
    % Cross-correlation using FFT-based convolution
    R = ifft(fft(received_signal_demod(n, :)) .* conj(fft(ref_signal)));
    range_profile(n, :) = abs(R); % Store magnitude
end

% Sum over antennas to improve SNR
range_profile_sum = sum(range_profile, 1);

% Find the peak corresponding to the target range
[~, idx_peak] = max(range_profile_sum);
tau_est = t(idx_peak); % Estimated time delay

% Estimate range
r_est = (tau_est * 3e8) / 2;

% Perform DoA Estimation at Estimated Range
signal_at_tau = received_signal_demod(:, idx_peak); % N x 1 vector

% Compute phase differences between adjacent antennas
delta_phi = angle(signal_at_tau(2:end) .* conj(signal_at_tau(1:end-1)));

% Average the phase differences
delta_phi_avg = mean(delta_phi);

% Estimate angle
theta_est_rad = asin((-delta_phi_avg * lambda) / (2 * pi * dx));
theta_est = rad2deg(theta_est_rad); 

% Output Results
fprintf('True Range: %.2f meters\n', r);
fprintf('Estimated Range: %.2f meters\n', r_est);
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est);

% Plot Range Profile
figure;
plot(t * 3e8 / 2, range_profile_sum, 'LineWidth', 1.5);
xlabel('Range (meters)');
ylabel('Amplitude');
title('Range Profile');
grid on;
xlim([r - 1, r + 1]);
hold on;
plot(r_est, range_profile_sum(idx_peak), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Range Profile', 'Estimated Range');

% Plot Phase Across Antennas at Estimated Range
figure;
plot(1:N, angle(signal_at_tau), 'o-', 'LineWidth', 1.5);
xlabel('Antenna Index');
ylabel('Phase (radians)');
title('Phase Across Antennas at Estimated Range');
grid on;


%% part 5-MIMO Array
clc; clear; close all;

% Parameters
c = 3e8;

f0 = 77e9;
lambda = c / f0;  
k = 2 * pi / lambda;  

% Antenna Array Parameters
NTx = 2;  
NRx = 58; 
N_virtual = NTx * NRx;  

% Virtual Array Spacing
dx = lambda/2; 

% Transmit Antenna Positions
Tx_positions = (0:NTx-1) * NRx * dx; % 1 x NTx

% Receive Antenna Positions
Rx_positions = (0:NRx-1) * dx; % 1 x NRx

% Target Parameters
theta_true = 28;  
theta_rad = deg2rad(theta_true); 
r = 78; % Target range 

% Time Delay Calculation
tau = 2 * r / c; % Round-trip time delay (seconds)

% Sampling Parameters
B = 1e9;  
Fs = 10 * B; % Sampling frequency (Hz) - oversampling by factor of 10
Ts = 1 / Fs; % Sampling interval 
t_max = tau + 5e-6; % Time window to capture the signal with margin
t = 0 : Ts : t_max; % Time vector

% Transmitted Signal Generation
g_tx = sinc(B * t) .* exp(1j * 2 * pi * f0 * t); % 1 x length(t)


% Initialize Received Signal Matrix (NTx x NRx x length(t))
received_signal = zeros(NTx, NRx, length(t));

% Expand Tx_positions and Rx_positions to matrices for pairwise computation
[Tx_matrix, Rx_matrix] = meshgrid(Tx_positions, Rx_positions);
Tx_matrix = Tx_matrix'; % NTx x NRx
Rx_matrix = Rx_matrix'; % NTx x NRx

% Calculate Phase Shifts for Each (Tx_i, Rx_j) Pair
phi_matrix = -k * (Tx_matrix + Rx_matrix) * sin(theta_rad); % NTx x NRx

% Generate Received Signals
for i = 1:NTx
    % Shift the transmitted signal by tau in samples
    shift_samples = round(tau / Ts); % Number of samples to shift
    if shift_samples < length(t)
        g_tx_shifted = [zeros(1, shift_samples), g_tx(1:end-shift_samples)];
    else
        g_tx_shifted = zeros(1, length(t));
    end
    
    for j = 1:NRx
        % Compute phase shift for Tx_i to Rx_j
        phi_ij = phi_matrix(i, j);
        
        % Apply phase shift
        g_rx = g_tx_shifted .* exp(1j * phi_ij);
        
        % Store the received signal at Rx_j from Tx_i
        received_signal(i, j, :) = g_rx;
    end
end


% Add Noise to Each Tx-Rx Pair Signal
SNR_dB = 30; 
signal_power = mean(abs(received_signal(:)).^2);
noise_power = signal_power / (10^(SNR_dB / 10));
noise = sqrt(noise_power/2) * (randn(NTx, NRx, length(t)) + 1j * randn(NTx, NRx, length(t)));
received_signal_noisy = received_signal + noise;

% Demodulate Each Tx-Rx Pair Signal to Baseband
carrier_conj = conj(exp(1j * 2 * pi * f0 * t)); % 1 x length(t)
received_signal_demod = received_signal_noisy .* reshape(carrier_conj, [1, 1, length(t)]); % NTx x NRx x length(t)

% Range Estimation via Matched Filtering

% Sum Over Transmit and Receive Antennas
received_signal_demod_sum = squeeze(sum(sum(received_signal_demod, 1), 2)); % length(t) x 1

% Transpose to make it a row vector
received_signal_demod_sum = received_signal_demod_sum.'; % 1 x length(t)

% Reference Signal for Matched Filtering (sinc pulse)
ref_signal = sinc(B * t); % 1 x length(t)

% Perform Matched Filtering
R = ifft(fft(received_signal_demod_sum) .* conj(fft(ref_signal)));

% Compute range profile
range_profile = abs(R);

% Find the Peak Corresponding to the Target Range
[~, idx_peak] = max(range_profile);
tau_est = t(idx_peak); % Estimated time delay

% Estimate Range
r_est = (tau_est * c) / 2; % Range estimation

% DoA Estimation Using Virtual Array

% Initialize Arrays for Virtual Array Positions and Signals
V_pos = zeros(N_virtual, 1);
V_signal = zeros(N_virtual, 1);

n = 1;
for i = 1:NTx
    for j = 1:NRx
        % Virtual Array Position
        V_pos(n) = Tx_positions(i) + Rx_positions(j);
        
        % Signal at Estimated Time Index
        V_signal(n) = received_signal_demod(i, j, idx_peak);
        
        n = n + 1;
    end
end

% Sort Virtual Array Positions and Signals
[V_pos_sorted, idx_sort] = sort(V_pos);
V_signal_sorted = V_signal(idx_sort);

% Verify Uniform Spacing
virtual_spacing = diff(V_pos_sorted);
mean_spacing = mean(virtual_spacing);
fprintf('Virtual array mean spacing: %.6e meters\n', mean_spacing);

% Adjust for Mean Spacing if Necessary
dx_virtual = mean_spacing; 

% Compute Phase Differences Between Adjacent Virtual Elements
delta_phi = angle(V_signal_sorted(2:end) .* conj(V_signal_sorted(1:end-1)));

% Unwrap Phase Differences
delta_phi_unwrapped = unwrap(delta_phi);

% Average the Phase Differences
delta_phi_avg = mean(delta_phi_unwrapped);

% Estimate Angle
theta_est_rad = asin((-delta_phi_avg * lambda) / (2 * pi * dx_virtual));
theta_est_deg = rad2deg(theta_est_rad); 

% Display Results

fprintf('True Range: %.2f meters\n', r);
fprintf('Estimated Range: %.2f meters\n', r_est);
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est_deg);



% Plot Range Profile
figure;
plot(t * c / 2, range_profile, 'LineWidth', 1.5);
xlabel('Range (meters)');
ylabel('Amplitude');
title('Range Profile');
grid on;
xlim([r - 1, r + 1]); 
hold on;
plot(r_est, range_profile(idx_peak), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Range Profile', 'Estimated Range');

% Plot Phase Across Virtual Array at Estimated Range
figure;
plot(V_pos_sorted, angle(V_signal_sorted), 'o-', 'LineWidth', 1.5);
xlabel('Virtual Array Position (meters)');
ylabel('Phase (radians)');
title('Phase Across Virtual Array at Estimated Range');
grid on;

%% part 6-Orthogonal Waveform
clc; clear; close all;

c = 3e8; 

f0 = 77e9;
lambda = c / f0; 
k = 2 * pi / lambda; % Wavenumber

% Antenna Array Parameters
NTx = 2; 
NRx = 4; 
N_virtual = NTx * NRx; 

% Virtual Array Spacing
dx = lambda / 2;  

% Transmit Antenna Positions
Tx_positions = (0:NTx-1) * NRx * dx; % 1 x NTx

% Receive Antenna Positions
Rx_positions = (0:NRx-1) * dx; % 1 x NRx

% Target Parameters
theta_true = 25; 
theta_rad = deg2rad(theta_true); 
r = 100; 

% Time Delay Calculation
tau = 2 * r / c; % Round-trip time delay (seconds)

% Sampling Parameters
B = 100e6; 
Fs = 2 * B;  
Ts = 1 / Fs; 
t_max = tau + 1e-6; 
t = 0 : Ts : t_max; 

% OFDM Parameters
N_subcarriers = NTx;  
subcarrier_spacing = B / N_subcarriers;  
subcarrier_frequencies = f0 + (-floor(N_subcarriers/2):floor((N_subcarriers-1)/2)) * subcarrier_spacing;

% Generate Orthogonal Waveforms (OFDM)

% Initialize transmitted signals matrix
g_tx = zeros(NTx, length(t)); % NTx x length(t)

for i = 1:NTx
    % Generate OFDM subcarrier signal for each transmitter
    g_tx(i, :) = exp(1j * 2 * pi * subcarrier_frequencies(i) * t);
end
 

% Initialize Received Signal Matrix (NTx x NRx x length(t))
received_signal = zeros(NTx, NRx, length(t));

% Expand Tx_positions and Rx_positions to matrices for pairwise computation
[Tx_matrix, Rx_matrix] = meshgrid(Tx_positions, Rx_positions);
Tx_matrix = Tx_matrix'; % NTx x NRx
Rx_matrix = Rx_matrix'; % NTx x NRx

% Calculate Phase Shifts for Each (Tx_i, Rx_j) Pair
phi_matrix = -k * (Tx_matrix + Rx_matrix) * sin(theta_rad); % NTx x NRx

% Generate Received Signals
for i = 1:NTx
    % Shift the transmitted signal by tau in samples
    shift_samples = round(tau / Ts); % Number of samples to shift
    if shift_samples < length(t)
        g_tx_shifted = [zeros(1, shift_samples), g_tx(i, 1:end-shift_samples)];
    else
        g_tx_shifted = zeros(1, length(t));
    end
    
    for j = 1:NRx
        % Compute phase shift for Tx_i to Rx_j
        phi_ij = phi_matrix(i, j);
        
        % Apply phase shift
        g_rx = g_tx_shifted .* exp(1j * phi_ij);
        
        % Store the received signal at Rx_j from Tx_i
        received_signal(i, j, :) = g_rx;
    end
end

% Add Noise to Each Tx-Rx Pair Signal

SNR_dB = 30;  
signal_power = mean(abs(received_signal(:)).^2);
noise_power = signal_power / (10^(SNR_dB / 10));
noise = sqrt(noise_power/2) * (randn(NTx, NRx, length(t)) + 1j * randn(NTx, NRx, length(t)));
received_signal_noisy = received_signal + noise;

% Matched Filtering  

% Initialize Matched Filter Outputs
matched_filter_output = zeros(NTx, NRx, length(t));

% Perform Matched Filtering for Each Transmit Signal
for i = 1:NTx
    % Reference Signal for Transmitter i
    ref_signal = conj(g_tx(i, :)); % 1 x length(t)
    
    for j = 1:NRx
        % Received Signal at Rx_j from Tx_i
        rx_signal = squeeze(received_signal_noisy(i, j, :)).'; % 1 x length(t)
        
        % Perform Matched Filtering
        R = ifft(fft(rx_signal) .* fft(ref_signal));
        
        % Store Matched Filter Output
        matched_filter_output(i, j, :) = R;
    end
end

% Range Estimation

% Sum over Receive Antennas for Each Transmitter
range_profile = zeros(NTx, length(t));

for i = 1:NTx
    % Sum over Receive Antennas
    range_profile(i, :) = sum(abs(matched_filter_output(i, :, :)), 2);
end

% Estimate Range for Each Transmitter
tau_est = zeros(NTx, 1);
r_est = zeros(NTx, 1);
idx_peak = zeros(NTx, 1);

for i = 1:NTx
    [~, idx_peak(i)] = max(range_profile(i, :));
    tau_est(i) = t(idx_peak(i));
    r_est(i) = (tau_est(i) * c) / 2;
end

% Use the average of estimated ranges
r_est_avg = mean(r_est);

% DoA Estimation Using Virtual Array

% Initialize Arrays for Virtual Array Positions and Signals
V_pos = zeros(N_virtual, 1);
V_signal = zeros(N_virtual, 1);

n = 1;
for i = 1:NTx
    for j = 1:NRx
        % Virtual Array Position
        V_pos(n) = Tx_positions(i) + Rx_positions(j);
        
        % Signal at Estimated Time Index
        V_signal(n) = matched_filter_output(i, j, idx_peak(i));
        
        n = n + 1;
    end
end

% Sort Virtual Array Positions and Signals
[V_pos_sorted, idx_sort] = sort(V_pos);
V_signal_sorted = V_signal(idx_sort);

% Compute Phase Differences Between Adjacent Virtual Elements
delta_phi = angle(V_signal_sorted(2:end) .* conj(V_signal_sorted(1:end-1)));

% Unwrap Phase Differences
delta_phi_unwrapped = unwrap(delta_phi);

% Average the Phase Differences
delta_phi_avg = mean(delta_phi_unwrapped);

% Estimate Angle
dx_virtual = mean(diff(V_pos_sorted)); % Use actual spacing
theta_est_rad = asin((-delta_phi_avg * lambda) / (2 * pi * dx_virtual));
theta_est_deg = rad2deg(theta_est_rad);  

% Display Results

fprintf('True Range: %.2f meters\n', r);
fprintf('Estimated Range: %.2f meters\n', r_est_avg);
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est_deg);

% Validate Orthogonality

% Cross-correlation between different transmit signals
cross_correlation = zeros(NTx, NTx);

for i = 1:NTx
    for j = 1:NTx
        if i ~= j
            cross_correlation(i, j) = max(abs(xcorr(g_tx(i, :), g_tx(j, :))));
        end
    end
end

fprintf('Maximum Cross-correlation between different transmit signals: %.2e\n', max(cross_correlation(:)));

% Plot Range Profiles for Each Transmitter

figure;
for i = 1:NTx
    plot(t * c / 2, range_profile(i, :), 'LineWidth', 1.5);
    hold on;
end
xlabel('Range (meters)');
ylabel('Amplitude');
title('Range Profiles for Each Transmitter');
legend('Tx1', 'Tx2');
grid on;

% Plot Phase Across Virtual Array at Estimated Range

figure;
plot(V_pos_sorted, angle(V_signal_sorted), 'o-', 'LineWidth', 1.5);
xlabel('Virtual Array Position (meters)');
ylabel('Phase (radians)');
title('Phase Across Virtual Array at Estimated Range');
grid on;
