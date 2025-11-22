% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

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

% Plot the sampled signals
figure('Name', 'Sampling x(t) at f_{Nyquist}');
plot(t, xt, '-', 'DisplayName', 'x(t)', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]); hold on;
stem(tSampleNqyuist, xSampleNyquist, 'DisplayName', 'x[n]', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (V)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('x(t) and x[n] sampled at f_{Nyquist}=(%.2f Hz)', fNyquist), 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name', 'Sampling x(t) at f_{s1}');
plot(t, xt, '-', 'DisplayName', 'x(t)', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]); hold on;
stem(tSample1, xSample1, 'DisplayName', 'x[n]', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (V)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('x(t) and x[n] sampled at f_{s1}=(%.2f Hz)', fS1), 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name', 'Sampling x(t) at f_{s2}');
plot(t, xt, '-', 'DisplayName', 'x(t)', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]); hold on;
stem(tSample2, xSample2, 'DisplayName', 'x[n]', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (V)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('x(t) and x[n] sampled at f_{s2}=(%.2f Hz)', fS2), 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Reconstruct the signals
xHatNyquist   = reconstruct(t, xSampleNyquist, fNyquist);
xHat1         = reconstruct(t, xSample1, fS1);
xHat2         = reconstruct(t, xSample2, fS2);

% Plot the reconstructed signals
figure('Name','Reconstruction From Samples at f_{Nyquist}');
plot(t, xt, '--', 'DisplayName', 'x(t) Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(t, xHatNyquist, '-', 'DisplayName', '\hat{x}(t) Reconstructed', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('x(t) reconstructed from x[n] sampled at f_{Nyquist}', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name','Reconstruction From Samples at f_{s1}');
plot(t, xt, '--', 'DisplayName', 'x(t) Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(t, xHat1, '-', 'DisplayName', '\hat{x}(t) Reconstructed', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('x(t) reconstructed from x[n] sampled at f_{s1}', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name','Reconstruction From Samples at f_{s2}');
plot(t, xt, '--', 'DisplayName', 'x(t) Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(t, xHat2, '-', 'DisplayName', '\hat{x}(t) Reconstructed', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (V)', 'FontSize', 12, 'FontWeight', 'bold');
title('x(t) reconstructed from x[n] sampled at f_{s2}', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Calculate the Mean Squared Error (MSE) for all three reconstructed signals
mseNyquist     = MSE(xt, xHatNyquist, t);
mseS1          = MSE(xt, xHat1, t);
mseS2          = MSE(xt, xHat2, t);

% Print the MSE values
fprintf('\nMean Squared Error from reconstruction at different frequencies:\n');
fprintf('* At f_{Nyquist}: %.3e\n', mseNyquist);
fprintf('* At f_{s1}: %.3e\n', mseS1);
fprintf('* At f_{s2}: %.3e\n', mseS2);
end