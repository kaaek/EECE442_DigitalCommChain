function [] = sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)

fprintf(repmat('=',1,70));
fprintf('\n   1.1 Generic Sampling is Running...\n');
fprintf(repmat('=',1,70));

[t, xt, f_max] = AMWave(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY);

fNyquist    = 2*f_max;
fS1         = 0.5*fNyquist;   % aliasing case
fS2         = 2*fNyquist;     % valid case

fprintf('\nMaximum frequency component in x(t): %.2f Hz\n', f_max);
fprintf('Sampling frequencies:\n');
fprintf('f_Nyquist = %.2f Hz\n', fNyquist);
fprintf('f_s1 = %.2f Hz\n', fS1);
fprintf('f_s2 = %.2f Hz\n', fS2);

% Sample to discrete x[n] (note that this is the analytical notation, but numerically, n is so fine we can still consider it continuous)

[tSampleNqyuist, xSampleNyquist]    = sample(t, xt, fNyquist);
[tSample1, xSample1]                = sample(t, xt, fS1);
[tSample2, xSample2]                = sample(t, xt, fS2);

% Reconstruct the signals
xHatNyquist   = reconstruct(t, xSampleNyquist, fNyquist);
xHat1         = reconstruct(t, xSample1, fS1);
xHat2         = reconstruct(t, xSample2, fS2);

% ===== Plot Sampling and Reconstruction at f_Nyquist =====
figure('Name', 'Sampling & Reconstruction at f_{Nyquist}', 'Position', [100 800 1200 500]);

% Subplot 1: Sampling at f_Nyquist
subplot(2, 1, 1);
plot(t, xt, '-', 'DisplayName', 'x(t)', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]); hold on;
stem(tSampleNqyuist, xSampleNyquist, 'DisplayName', 'Samples', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title(sprintf('Sampling: f_s = %.0f Hz', fNyquist), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% Subplot 2: Reconstruction at f_Nyquist
subplot(2, 1, 2);
plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(t, xHatNyquist, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title(sprintf('Reconstruction: f_s = %.0f Hz', fNyquist), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% ===== Plot Sampling and Reconstruction at f_s1 =====
figure('Name', 'Sampling & Reconstruction at f_{s1}', 'Position', [100 400 1200 500]);

% Subplot 1: Sampling at f_s1
subplot(2, 1, 1);
plot(t, xt, '-', 'DisplayName', 'x(t)', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]); hold on;
stem(tSample1, xSample1, 'DisplayName', 'Samples', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title(sprintf('Sampling: f_s = %.0f Hz (Undersampled)', fS1), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% Subplot 2: Reconstruction at f_s1
subplot(2, 1, 2);
plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(t, xHat1, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title(sprintf('Reconstruction: f_s = %.0f Hz (Aliased)', fS1), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% ===== Plot Sampling and Reconstruction at f_s2 =====
figure('Name', 'Sampling & Reconstruction at f_{s2}', 'Position', [100 0 1200 500]);

% Subplot 1: Sampling at f_s2
subplot(2, 1, 1);
plot(t, xt, '-', 'DisplayName', 'x(t)', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]); hold on;
stem(tSample2, xSample2, 'DisplayName', 'Samples', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title(sprintf('Sampling: f_s = %.0f Hz (Oversampled)', fS2), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% Subplot 2: Reconstruction at f_s2
subplot(2, 1, 2);
plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(t, xHat2, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title(sprintf('Reconstruction: f_s = %.0f Hz (Perfect)', fS2), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 9);

% Calculate the Mean Squared Error (MSE) for all three reconstructed signals
mseNyquist     = MSE(xt, xHatNyquist, t);
mseS1          = MSE(xt, xHat1, t);
mseS2          = MSE(xt, xHat2, t);

% Print the MSE values
fprintf('\nMean Squared Error from reconstruction at different frequencies:\n');
fprintf('* At f_{Nyquist}: %.3e\n', mseNyquist);
fprintf('* At f_{s1}: %.3e\n', mseS1);
fprintf('* At f_{s2}: %.3e\n', mseS2);

figs = findall(0, 'Type', 'figure');
set(figs, 'Units', 'pixels', 'Position', [100 100 800 600]);
end