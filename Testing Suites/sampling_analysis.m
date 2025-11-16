% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [] = sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)

fprintf('========================================\n');
fprintf('   1.1 Generic Sampling is Running...\n');
fprintf('========================================\n');

[t, xt, f_max] = AMWave(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY);

fNyquist    = 2*f_max;
fS1         = 0.5*fNyquist;   % aliasing case
fS2         = 2*fNyquist;     % valid case

% Sample to discrete x[n] (note that this is the analytical notation, but numerically, n is so fine we can still consider it continuous)

[tSampleNqyuist, xSampleNyquist]    = sample(t, xt, fNyquist);
[tSample1, xSample1]                = sample(t, xt, fS1);
[tSample2, xSample2]                = sample(t, xt, fS2);

% Plot the sampled signals
figure('Name', 'Sampler');

subplot(3, 1, 1);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
stem(tSampleNqyuist, xSampleNyquist, 'b', 'DisplayName', 'x[n] sampled at f_{Nyquist}', 'Marker', 'o');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title(sprintf('Comparison of x(t) and x[n] sampled at f_{Nyquist} (%.2f Hz)', fNyquist));
legend show;
grid on;

subplot(3, 1, 2);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
stem(tSample1, xSample1, 'g', 'DisplayName', 'x[n] sampled at f_{s1}', 'Marker', 'x');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title(sprintf('Comparison of x(t) and x[n] sampled at f_{s1} (%.2f Hz)', fS1));
legend show;
grid on;

subplot(3, 1, 3);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
stem(tSample2, xSample2, 'm', 'DisplayName', 'x[n] sampled at f_{s2}', 'Marker', 's');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title(sprintf('Comparison of x(t) and x[n] sampled at f_{s2} (%.2f Hz)', fS2));
legend show;
grid on;

% Reconstruct the signals
xHatNyquist   = reconstruct(t, xSampleNyquist, fNyquist);
xHat1         = reconstruct(t, xSample1, fS1);
xHat2         = reconstruct(t, xSample2, fS2);

% Plot the reconstructed signals
figure('Name','Reconstruction From Samples');

subplot(3, 1, 1);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
plot(t, xHatNyquist, 'b', 'DisplayName', 'Reconstructed x[n] at f_{Nyquist}', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Reconstruction of x[n] sampled at f_{Nyquist}');
legend show;
grid on;

subplot(3, 1, 2);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
plot(t, xHat1, 'g', 'DisplayName', 'Reconstructed x[n] at f_{s1}', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Reconstruction of x[n] sampled at f_{s1}');
legend show;
grid on;

subplot(3, 1, 3);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
plot(t, xHat2, 'm', 'DisplayName', 'Reconstructed x[n] at f_{s2}', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Reconstruction of x[n] sampled at f_{s2}');
legend show;
grid on;

% Calculate the Mean Squared Error (MSE) for all three reconstructed signals
mseNyquist     = MSE(xt, xHatNyquist, t);
mseS1          = MSE(xt, xHat1, t);
mseS2          = MSE(xt, xHat2, t);

% Print the MSE values
fprintf('Mean Squared Error for reconstruction:\n');
fprintf('* At f_{Nyquist}: %.3e\n', mseNyquist);
fprintf('* At f_{s1}: %.3e\n', mseS1);
fprintf('* At f_{s2}: %.3e\n', mseS2);
end