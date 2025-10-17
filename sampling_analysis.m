% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function sampling_analysis(SIGNALDURATION, CARRIERFREQUENCY)

fprintf('========================================\n');
fprintf('   1.1 Generic Sampling is Running...\n');
fprintf('========================================\n');

[t, xt, f_max] = exampleSpeechWave(SIGNALDURATION, CARRIERFREQUENCY);
f_Nyquist = 2*f_max;
f_s1 = 0.5*f_Nyquist;   % aliasing case
f_s2 = 2*f_Nyquist;     % valid case

% Sample to discrete x[n] (note that this is the analytical notation, but
% numerically, n is so fine we can still consider it continuous)

[t_sample_nqyuist, x_sample_nyquist] = sample(t, xt, f_Nyquist);
[t_sample1, x_sample1] = sample(t, xt, f_s1);
[t_sample2, x_sample2] = sample(t, xt, f_s2);

% Plot the sampled signals
figure('Name', 'Sampler');

subplot(3, 1, 1);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
stem(t_sample_nqyuist, x_sample_nyquist, 'b', 'DisplayName', 'x[n] sampled at f_{Nyquist}', 'Marker', 'o');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title(sprintf('Comparison of x(t) and x[n] sampled at f_{Nyquist} (%.2f Hz)', f_Nyquist));
legend show;
grid on;

subplot(3, 1, 2);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
stem(t_sample1, x_sample1, 'g', 'DisplayName', 'x[n] sampled at f_{s1}', 'Marker', 'x');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title(sprintf('Comparison of x(t) and x[n] sampled at f_{s1} (%.2f Hz)', f_s1));
legend show;
grid on;

subplot(3, 1, 3);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
stem(t_sample2, x_sample2, 'm', 'DisplayName', 'x[n] sampled at f_{s2}', 'Marker', 's');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title(sprintf('Comparison of x(t) and x[n] sampled at f_{s2} (%.2f Hz)', f_s2));
legend show;
grid on;

% Reconstruct the signals
x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2  = reconstruct(t, x_sample2, f_s2);

% Plot the reconstructed signals
figure('Name','Reconstruction From Samples');

subplot(3, 1, 1);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
plot(t, x_hat_nyquist, 'b', 'DisplayName', 'Reconstructed x[n] at f_{Nyquist}', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Reconstruction of x[n] sampled at f_{Nyquist}');
legend show;
grid on;

subplot(3, 1, 2);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
plot(t, x_hat_1, 'g', 'DisplayName', 'Reconstructed x[n] at f_{s1}', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Reconstruction of x[n] sampled at f_{s1}');
legend show;
grid on;

subplot(3, 1, 3);
plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]); hold on;
plot(t, x_hat_2, 'm', 'DisplayName', 'Reconstructed x[n] at f_{s2}', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Reconstruction of x[n] sampled at f_{s2}');
legend show;
grid on;

% Calculate the Mean Squared Error (MSE) for all three reconstructed signals
mse_nyquist = MSE(t, xt, x_hat_nyquist);
mse_s1 = MSE(t, xt, x_hat_1);
mse_s2 = MSE(t, xt, x_hat_2);

% Print the MSE values
fprintf('Mean Squared Error for reconstruction:\n');
fprintf('* At f_{Nyquist}: %.3e\n', mse_nyquist);
fprintf('* At f_{s1}: %.3e\n', mse_s1);
fprintf('* At f_{s2}: %.3e\n', mse_s2);

% TO-DO: 
% * Commentary on which reconstruction better matches and why (reference sampling theorem).
% * Show aliasing visually.
% * Discuss frequency-domain folding due to undersampling.
end