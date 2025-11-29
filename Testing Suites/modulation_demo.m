% EXAMPLE 1: BPSK Modulation and Demodulation
fprintf('EXAMPLE 1: BPSK Modulation/Demodulation\n');
fprintf(repmat('-', 1, 70));

num_bits_bpsk = 16; % Generate a random bitstream for BPSK
bitstream_bpsk = randi([0, 1], 1, num_bits_bpsk);
fprintf('\nInput bitstream (BPSK): %s\n', sprintf('%d', bitstream_bpsk));

% BPSK Modulation (apply to each bit)
% Using A = sqrt(2) for unit energy symbols
modulated_bpsk = zeros(1, num_bits_bpsk);
A_bpsk = sqrt(2);  % Amplitude for unit energy normalization
for i = 1:num_bits_bpsk
    modulated_bpsk(i) = bpsk_mod(bitstream_bpsk(i), A_bpsk);  % normalized amplitude
end
fprintf('Modulated BPSK symbols: %s\n', sprintf('%+.1f ', real(modulated_bpsk)));

% Create transmit signal (sinusoidal carrier)
fs = 1000;  % Sampling frequency
fc = 100;   % Carrier frequency
t_bpsk = (0:length(modulated_bpsk)-1) / fs;
carrier_bpsk = exp(1j * 2 * pi * fc * t_bpsk);
tx_signal_bpsk = modulated_bpsk .* carrier_bpsk;

% BPSK Demodulation (with noise-free channel)
demod_bpsk = zeros(1, num_bits_bpsk);
for i = 1:num_bits_bpsk
    demod_bpsk(i) = bpsk_demod(modulated_bpsk(i));
end
fprintf('Demodulated bitstream: %s\n', sprintf('%d', demod_bpsk));
fprintf('Bit Error Rate (BER): %.2f%%\n', (sum(bitstream_bpsk ~= demod_bpsk) / num_bits_bpsk) * 100);

% Visualize BPSK
figure('Name', 'BPSK Modulation Example');

% Constellation
subplot(2, 2, 1);
scatter(real(modulated_bpsk), imag(modulated_bpsk), 100, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
xlabel('In-Phase (I)', 'FontSize', 11);
ylabel('Quadrature (Q)', 'FontSize', 11);
title('BPSK Constellation', 'FontSize', 12);
grid on;
axis([-1.5 1.5 -1.5 1.5]);
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Time-domain modulated signal
subplot(2, 2, 2);
plot(t_bpsk(1:min(100, length(t_bpsk))), real(tx_signal_bpsk(1:min(100, length(tx_signal_bpsk)))), 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Amplitude', 'FontSize', 11);
title('BPSK Modulated Signal (Real Part)', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Bitstream representation
subplot(2, 2, 3);
stem(0:length(bitstream_bpsk)-1, bitstream_bpsk, 'b', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', [0 0.45 0.74]);
xlabel('Bit Index', 'FontSize', 11);
ylabel('Bit Value', 'FontSize', 11);
title('Input Bitstream', 'FontSize', 12);
ylim([-0.5 1.5]);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Constellation mapping
subplot(2, 2, 4);
scatter(1, 0, 200, [0 0.45 0.74], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
hold on;
scatter(-1, 0, 200, [0.85 0.33 0.1], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
text(1, 0.25, '1', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(-1, 0.25, '0', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('In-Phase (I)', 'FontSize', 11);
ylabel('Quadrature (Q)', 'FontSize', 11);
title('BPSK Symbol Mapping', 'FontSize', 12);
axis([-1.5 1.5 -0.5 0.5]);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

fprintf('\n');

% EXAMPLE 2: QPSK Modulation and Demodulation
fprintf('EXAMPLE 2: QPSK Modulation/Demodulation\n');
fprintf(repmat('-', 1, 70));

% Generate a random bitstream for QPSK (must be even length)
num_bits_qpsk = 16;
bitstream_qpsk = randi([0, 1], 1, num_bits_qpsk);
fprintf('\nInput bitstream (QPSK): %s\n', sprintf('%d', bitstream_qpsk));

% QPSK Modulation (pairs of bits at a time)
% Using A = sqrt(2) for unit energy symbols (standard in digital communications)
num_symbols_qpsk = num_bits_qpsk / 2;
modulated_qpsk = zeros(1, num_symbols_qpsk);
symbol_strings = cell(1, num_symbols_qpsk);
A_qpsk = sqrt(2);  % Amplitude for unit energy normalization
for i = 1:num_symbols_qpsk
    bit_pair = sprintf('%d%d', bitstream_qpsk(2*i-1), bitstream_qpsk(2*i));
    symbol_strings{i} = bit_pair;
    modulated_qpsk(i) = qpsk_mod(bit_pair, A_qpsk);  % Normalized amplitude
end
fprintf('Modulated QPSK symbols (Gray-coded): %s\n', ...
    sprintf('%.2f%+.2fj ', real(modulated_qpsk), imag(modulated_qpsk)));

% Create transmit signal (sinusoidal carrier)
t_qpsk = (0:length(modulated_qpsk)-1) / fs;
carrier_qpsk = exp(1j * 2 * pi * fc * t_qpsk);
tx_signal_qpsk = modulated_qpsk .* carrier_qpsk;

% QPSK Demodulation (with noise-free channel) - check symbol-level fidelity
fprintf('Demodulated QPSK symbols (should match input):\n');
fprintf('Symbol # | Modulated Symbol | I | Q | Decoded\n');
fprintf(repmat('-', 1, 60));
fprintf('\n');
for i = 1:num_symbols_qpsk
    rx_i = real(modulated_qpsk(i));
    rx_q = imag(modulated_qpsk(i));
    
    % Decode based on quadrants
    if rx_i >= 0 && rx_q >= 0
        decoded = '11';
    elseif rx_i < 0 && rx_q >= 0
        decoded = '10';
    elseif rx_i < 0 && rx_q < 0
        decoded = '00';
    else  % rx_i >= 0 && rx_q < 0
        decoded = '01';
    end
    
    fprintf('   %d    | %+.2f%+.2fj      | %+.2f | %+.2f | %s (expected: %s)\n', ...
        i, rx_i, rx_q, rx_i, rx_q, decoded, symbol_strings{i});
end
fprintf('\n');

% Visualize QPSK
figure('Name', 'QPSK Modulation Example');

% Constellation
subplot(2, 2, 1);
scatter(real(modulated_qpsk), imag(modulated_qpsk), 100, [0 0.45 0.74], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
xlabel('In-Phase (I)', 'FontSize', 11);
ylabel('Quadrature (Q)', 'FontSize', 11);
title('QPSK Constellation', 'FontSize', 12);
grid on;
axis([-1.5 1.5 -1.5 1.5]);
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Time-domain modulated signal
subplot(2, 2, 2);
plot(t_qpsk(1:min(100, length(t_qpsk))), real(tx_signal_qpsk(1:min(100, length(tx_signal_qpsk)))), ...
     'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Amplitude', 'FontSize', 11);
title('QPSK Modulated Signal (Real Part)', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Bitstream representation
subplot(2, 2, 3);
stem(0:length(bitstream_qpsk)-1, bitstream_qpsk, 'b', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', [0 0.45 0.74]);
xlabel('Bit Index', 'FontSize', 11);
ylabel('Bit Value', 'FontSize', 11);
title('Input Bitstream', 'FontSize', 12);
ylim([-0.5 1.5]);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Constellation mapping (Counter-clockwise from '11')
subplot(2, 2, 4);
angles = [pi/4, 3*pi/4, 5*pi/4, 7*pi/4];
labels = {'11', '10', '00', '01'};  % Counter-clockwise order
colors = [0 0.45 0.74; 0.85 0.33 0.1; 0.47 0.67 0.19; 0.93 0.69 0.13];
for k = 1:4
    x = cos(angles(k))/sqrt(2);
    y = sin(angles(k))/sqrt(2);
    scatter(x, y, 200, colors(k, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    hold on;
    text(x*1.4, y*1.4, labels{k}, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');
end
xlabel('In-Phase (I)', 'FontSize', 11);
ylabel('Quadrature (Q)', 'FontSize', 11);
title('QPSK Symbol Mapping (Counter-clockwise)', 'FontSize', 12);
axis([-1.2 1.2 -1.2 1.2]);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

fprintf('\n');

% EXAMPLE 3: BPSK with AWGN Noise
fprintf('EXAMPLE 3: BPSK with AWGN Noise\n');
fprintf(repmat('-', 1, 70));

num_bits_noisy = 64;
bitstream_noisy = randi([0, 1], 1, num_bits_noisy);
modulated_noisy = zeros(1, num_bits_noisy);
A_bpsk_noisy = sqrt(2);
for i = 1:num_bits_noisy
    modulated_noisy(i) = bpsk_mod(bitstream_noisy(i), A_bpsk_noisy);
end

% Add AWGN noise at different SNR levels
snr_db = [0, 5, 10];

figure('Name', 'BPSK with AWGN Noise');

for snr_idx = 1:length(snr_db)
    snr = snr_db(snr_idx);
    signal_power = mean(abs(modulated_noisy).^2);
    noise_power = signal_power / (10^(snr/10));
    noise = sqrt(noise_power) * (randn(size(modulated_noisy)) + 1j * randn(size(modulated_noisy))) / sqrt(2);
    received_signal = modulated_noisy + noise;
    
    % Demodulate
    demod_noisy = zeros(1, num_bits_noisy);
    for i = 1:num_bits_noisy
        demod_noisy(i) = bpsk_demod(received_signal(i));
    end
    ber = sum(bitstream_noisy ~= demod_noisy) / num_bits_noisy;
    
    fprintf('\nSNR = %2d dB: BER = %.4f (%d/%d errors)\n', snr, ber, ...
            sum(bitstream_noisy ~= demod_noisy), num_bits_noisy);
    
    % Plot constellation
    subplot(1, 3, snr_idx);
    scatter(real(received_signal), imag(received_signal), 50, [0 0.45 0.74], ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
    hold on;
    scatter([1, -1], [0, 0], 150, [0.85 0.33 0.1], 'x', 'LineWidth', 2);
    xlabel('In-Phase (I)', 'FontSize', 11);
    ylabel('Quadrature (Q)', 'FontSize', 11);
    title(sprintf('SNR = %d dB (BER = %.4f)', snr, ber), 'FontSize', 12);
    axis([-2 2 -2 2]);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

fprintf('\n');

% EXAMPLE 4: QPSK with AWGN Noise
fprintf('EXAMPLE 4: QPSK with AWGN Noise\n');
fprintf(repmat('-', 1, 70));

num_bits_qpsk_noisy = 64;
bitstream_qpsk_noisy = randi([0, 1], 1, num_bits_qpsk_noisy);
num_symbols_qpsk_noisy = num_bits_qpsk_noisy / 2;
modulated_qpsk_noisy = zeros(1, num_symbols_qpsk_noisy);
A_qpsk_noisy = sqrt(2);
for i = 1:num_symbols_qpsk_noisy
    bit_pair = sprintf('%d%d', bitstream_qpsk_noisy(2*i-1), bitstream_qpsk_noisy(2*i));
    modulated_qpsk_noisy(i) = qpsk_mod(bit_pair, A_qpsk_noisy);
end

figure('Name', 'QPSK with AWGN Noise');

for snr_idx = 1:length(snr_db)
    snr = snr_db(snr_idx);
    signal_power = mean(abs(modulated_qpsk_noisy).^2);
    noise_power = signal_power / (10^(snr/10));
    noise = sqrt(noise_power) * (randn(size(modulated_qpsk_noisy)) + 1j * randn(size(modulated_qpsk_noisy))) / sqrt(2);
    received_signal_qpsk = modulated_qpsk_noisy + noise;
    
    % Demodulate - hard-decision decoder for each symbol
    % Use simple quadrant decision (standard QPSK receiver)
    num_errors = 0;
    for i = 1:num_symbols_qpsk_noisy
        bit_pair = sprintf('%d%d', bitstream_qpsk_noisy(2*i-1), bitstream_qpsk_noisy(2*i));
        expected_sym = qpsk_mod(bit_pair, 1);
        
        % Hard-decision decoding: check quadrants
        rx_i = real(received_signal_qpsk(i));
        rx_q = imag(received_signal_qpsk(i));
        
        % Decode based on quadrants
        if rx_i >= 0 && rx_q >= 0
            decoded = '11';  % corresponds to angle pi/4
        elseif rx_i < 0 && rx_q >= 0
            decoded = '10';  % corresponds to angle 3*pi/4
        elseif rx_i < 0 && rx_q < 0
            decoded = '00';  % corresponds to angle 5*pi/4
        else  % rx_i >= 0 && rx_q < 0
            decoded = '01';  % corresponds to angle 7*pi/4
        end
        
        % Compare with expected bits
        if ~strcmp(decoded, bit_pair)
            num_errors = num_errors + sum(str2num(decoded)' ~= str2num(bit_pair)');
        end
    end
    ber = num_errors / num_bits_qpsk_noisy;
    
    fprintf('\nSNR = %2d dB: BER = %.4f (%d/%d bit errors)\n', snr, ber, num_errors, num_bits_qpsk_noisy);
    
    % Plot constellation
    subplot(1, 3, snr_idx);
    scatter(real(received_signal_qpsk), imag(received_signal_qpsk), 50, [0 0.45 0.74], ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
    hold on;
    % Ideal QPSK points
    angles = [pi/4, 3*pi/4, 5*pi/4, 7*pi/4];
    for k = 1:4
        scatter(cos(angles(k))/sqrt(2), sin(angles(k))/sqrt(2), 150, [0.85 0.33 0.1], 'x', 'LineWidth', 2);
    end
    xlabel('In-Phase (I)', 'FontSize', 11);
    ylabel('Quadrature (Q)', 'FontSize', 11);
    title(sprintf('SNR = %d dB (BER = %.4f)', snr, ber), 'FontSize', 12);
    axis([-1.2 1.2 -1.2 1.2]);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\nDemonstration complete!\n');
fprintf(repmat('=', 1, 70));
