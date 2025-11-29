% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% ----------------------------------------------------------------------

function quantization_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)

fprintf(repmat('=',1,70));
fprintf('\n   2.1 Quantization Analysis running...\n');
fprintf(repmat('=',1,70));

% The below sampling is quiet since this is not the suite for sampling
[t, xt, f_max]  = AMWave(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY);
f_Nyquist       = 2*f_max;
f_s1            = 0.5*f_Nyquist;   % aliasing case
f_s2            = 2*f_Nyquist;     % valid case
% Sample to discrete x[n] (note that this is the analytical notation, but
% numerically, n is so fine we can still consider it continuous)
[t_sample_nqyuist, x_sample_nyquist]    = sample(t, xt, f_Nyquist);
[t_sample1, x_sample1]                  = sample(t, xt, f_s1);
[t_sample2, x_sample2]                  = sample(t, xt, f_s2);

% Quantizer parameters
M1 = 4;
M2 = 16;
M3 = 32;

% Two-level quantizer
xqn_2lvl = twoLvlQuan(x_sample_nyquist);
xq1_2lvl = twoLvlQuan(x_sample1);
xq2_2lvl = twoLvlQuan(x_sample2);

% Uniform quantization
[xqn_u1, ~, ~, ~] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M1);
[xq1_u1, ~, ~, ~] = uniformQuan(t_sample1, x_sample1, M1);
[xq2_u1, ~, ~, ~] = uniformQuan(t_sample2, x_sample2, M1);

[xqn_u2, ~, ~, ~] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M2);
[xq1_u2, ~, ~, ~] = uniformQuan(t_sample1, x_sample1, M2);
[xq2_u2, ~, ~, ~] = uniformQuan(t_sample2, x_sample2, M2);

[xqn_u3, ~, ~, ~] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M3);
[xq1_u3, ~, ~, ~] = uniformQuan(t_sample1, x_sample1, M3);
[xq2_u3, ~, ~, ~] = uniformQuan(t_sample2, x_sample2, M3);

% Lloyd-Max quantization
[thr_n_l1, lvl_n_l1, xqn_l1, ~] = lloydMax(x_sample_nyquist, M1, 0.1);
[thr_1_l1, lvl_1_l1, xq1_l1, ~] = lloydMax(x_sample1, M1, 0.1);
[thr_2_l1, lvl_2_l1, xq2_l1, ~] = lloydMax(x_sample2, M1, 0.1);

[thr_n_l2, lvl_n_l2, xqn_l2, ~] = lloydMax(x_sample_nyquist, M2, 0.1);
[thr_1_l2, lvl_1_l2, xq1_l2, ~] = lloydMax(x_sample1, M2, 0.1);
[thr_2_l2, lvl_2_l2, xq2_l2, ~] = lloydMax(x_sample2, M2, 0.1);

[thr_n_l3, lvl_n_l3, xqn_l3, ~] = lloydMax(x_sample_nyquist, M3, 0.1);
[thr_1_l3, lvl_1_l3, xq1_l3, ~] = lloydMax(x_sample1, M3, 0.1);
[thr_2_l3, lvl_2_l3, xq2_l3, ~] = lloydMax(x_sample2, M3, 0.1);

% Reconstruct all signals
x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2 = reconstruct(t, x_sample2, f_s2);

% Create figures organized by M value: for each M, create 3 figures (one per frequency)

% Two-level quantizer: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_2lvl, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Two-Level Quantizer', 2);
plotSingleQuantizer(t_sample1, x_sample1, xq1_2lvl, t, xt, x_hat_1, f_s1, 'Undersampled', 'Two-Level Quantizer', 2);
plotSingleQuantizer(t_sample2, x_sample2, xq2_2lvl, t, xt, x_hat_2, f_s2, 'Oversampled', 'Two-Level Quantizer', 2);

% Uniform M=4: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_u1, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Uniform Quantizer', 4);
plotSingleQuantizer(t_sample1, x_sample1, xq1_u1, t, xt, x_hat_1, f_s1, 'Undersampled', 'Uniform Quantizer', 4);
plotSingleQuantizer(t_sample2, x_sample2, xq2_u1, t, xt, x_hat_2, f_s2, 'Oversampled', 'Uniform Quantizer', 4);

% Uniform M=16: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_u2, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Uniform Quantizer', 16);
plotSingleQuantizer(t_sample1, x_sample1, xq1_u2, t, xt, x_hat_1, f_s1, 'Undersampled', 'Uniform Quantizer', 16);
plotSingleQuantizer(t_sample2, x_sample2, xq2_u2, t, xt, x_hat_2, f_s2, 'Oversampled', 'Uniform Quantizer', 16);

% Uniform M=32: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_u3, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Uniform Quantizer', 32);
plotSingleQuantizer(t_sample1, x_sample1, xq1_u3, t, xt, x_hat_1, f_s1, 'Undersampled', 'Uniform Quantizer', 32);
plotSingleQuantizer(t_sample2, x_sample2, xq2_u3, t, xt, x_hat_2, f_s2, 'Oversampled', 'Uniform Quantizer', 32);

% Lloyd-Max M=4: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_l1, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Lloyd-Max Quantizer', 4);
plotSingleQuantizer(t_sample1, x_sample1, xq1_l1, t, xt, x_hat_1, f_s1, 'Undersampled', 'Lloyd-Max Quantizer', 4);
plotSingleQuantizer(t_sample2, x_sample2, xq2_l1, t, xt, x_hat_2, f_s2, 'Oversampled', 'Lloyd-Max Quantizer', 4);

% Lloyd-Max M=16: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_l2, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Lloyd-Max Quantizer', 16);
plotSingleQuantizer(t_sample1, x_sample1, xq1_l2, t, xt, x_hat_1, f_s1, 'Undersampled', 'Lloyd-Max Quantizer', 16);
plotSingleQuantizer(t_sample2, x_sample2, xq2_l2, t, xt, x_hat_2, f_s2, 'Oversampled', 'Lloyd-Max Quantizer', 16);

% Lloyd-Max M=32: 3 figures (Nyquist, Undersampled, Oversampled)
plotSingleQuantizer(t_sample_nqyuist, x_sample_nyquist, xqn_l3, t, xt, x_hat_nyquist, f_Nyquist, 'Nyquist', 'Lloyd-Max Quantizer', 32);
plotSingleQuantizer(t_sample1, x_sample1, xq1_l3, t, xt, x_hat_1, f_s1, 'Undersampled', 'Lloyd-Max Quantizer', 32);
plotSingleQuantizer(t_sample2, x_sample2, xq2_l3, t, xt, x_hat_2, f_s2, 'Oversampled', 'Lloyd-Max Quantizer', 32);

% -------------------------------------------------------------------

% * Describe the visible stair-step distortion.
% * Comment on binary-level use cases (hard decisions, low SNR scenarios)

M = 1:1:256;
MSE = zeros(1, length(M));
for i = 1:length(M)
    % Perform uniform quantization for each M value
    [~, ~, ~, MSE_i] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M(i));
    MSE(i) = MSE_i;
end
figure('Name', 'Uniform Quantizer: MSE vs M');
plot(M, MSE, 'LineWidth', 2, 'Color', [0 0.45 0.74]);
xlabel('Number of Quantization Levels (M)', 'FontSize', 11);
ylabel('Mean Squared Error', 'FontSize', 11);
title('MSE vs Quantization Levels', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% PDF visualization for Lloyd-Max quantizers
plotPdf(x_sample_nyquist, thr_n_l1, lvl_n_l1, M1, 'Nyquist');
plotPdf(x_sample1, thr_1_l1, lvl_1_l1, M1, 'Undersampled');
plotPdf(x_sample2, thr_2_l1, lvl_2_l1, M1, 'Oversampled');

plotPdf(x_sample_nyquist, thr_n_l2, lvl_n_l2, M2, 'Nyquist');
plotPdf(x_sample1, thr_1_l2, lvl_1_l2, M2, 'Undersampled');
plotPdf(x_sample2, thr_2_l2, lvl_2_l2, M2, 'Oversampled');

plotPdf(x_sample_nyquist, thr_n_l3, lvl_n_l3, M3, 'Nyquist');
plotPdf(x_sample1, thr_1_l3, lvl_1_l3, M3, 'Undersampled');
plotPdf(x_sample2, thr_2_l3, lvl_2_l3, M3, 'Oversampled');


end

function plotSingleQuantizer(t_samp, x_samp, xq, t_rec, xt_rec, x_hat, f_s, freq_label, quantizer_name, M)
    % PLOTSINGLEQUANTIZER Creates a figure with quantized and reconstructed signals stacked vertically
    %
    % Inputs:
    %   t_samp, x_samp - Sampled time and signal
    %   xq - Quantized signal
    %   t_rec, xt_rec - Reconstruction time and original signal
    %   x_hat - Reconstructed signal
    %   f_s - Sampling frequency
    %   freq_label - Frequency label (e.g., 'Nyquist', 'Undersampled', 'Oversampled')
    %   quantizer_name - Name of quantizer
    %   M - Quantization levels
    
    figure('Name', sprintf('%s M=%d - %s (f_s=%.2f Hz)', quantizer_name, M, freq_label, f_s));
    
    % Quantized signal (top subplot)
    subplot(2, 1, 1);
    plot(t_samp, x_samp, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
    hold on;
    stem(t_samp, xq, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, ...
         'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Amplitude', 'FontSize', 10);
    title('Quantized Signal', 'FontSize', 11);
    legend('show', 'Location', 'best', 'FontSize', 8);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
    
    % Reconstructed signal (bottom subplot)
    subplot(2, 1, 2);
    plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
    hold on;
    plot(t_rec, x_hat, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Amplitude', 'FontSize', 10);
    title('Reconstructed Signal', 'FontSize', 11);
    legend('show', 'Location', 'best', 'FontSize', 8);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
end

function plotPdf(x_samples, thr, lvl, M, freq_label)
    % PLOTPDF Plots the empirical probability density function (PDF) of the input samples
    %
    % This function generates a plot of the empirical PDF of the input 
    % samples along with the quantization thresholds and representation 
    % levels. It visualizes how the quantization process affects the 
    % distribution of the input data.
    %
    % Inputs:
    %   x_samples - A vector of input samples for which the PDF is to be plotted.
    %   thr       - A vector of thresholds used for quantization.
    %   lvl       - A vector of representation levels after quantization.
    %   M         - Number of quantization levels.
    %   freq_label - Frequency case label ('Nyquist', 'Undersampled', 'Oversampled')
    %
    % Example:
    %   plotPdf(x_samples, thr, lvl, M, 'Nyquist');
    
    % Get the empirical PDF from the sampled sequence x_samples
    x_samples = real(x_samples);  % Disregard imaginary part
    nbins = min(max(80, round(length(x_samples) / 10)), 200);
    [pdf, edges] = histcounts(x_samples, nbins, 'Normalization', 'pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    figure('Name', sprintf('Lloyd-Max PDF M=%d - %s', M, freq_label));
    plot(centers, pdf, 'LineWidth', 2, 'Color', [0 0.45 0.74], 'DisplayName', 'PDF');
    hold on;
    for k = 1:length(thr)
        plot([thr(k) thr(k)], ylim, '--', 'LineWidth', 1.5, 'Color', [0.93 0.69 0.13]);
    end
    stem(lvl, interp1(centers, pdf, lvl, 'linear', 'extrap'), 'filled', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1], 'DisplayName', 'Levels');
    
    % Add threshold line to legend (just the first one for clarity)
    plot(thr(1), 0, '--', 'LineWidth', 1.5, 'Color', [0.93 0.69 0.13], 'DisplayName', 'Thresholds');
    
    legend('show', 'Location', 'best', 'FontSize', 8);
    title(sprintf('Lloyd-Max Quantization PDF (M=%d, %s)', M, freq_label), 'FontSize', 11);
    xlabel('Amplitude', 'FontSize', 10);
    ylabel('Probability Density', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
end

% function plotQuantizerByMValue(t_n, x_n, xq_n, t_1, x_1, xq_1, t_2, x_2, xq_2, ...
%                                 t_rec, xt_rec, x_hat_n, x_hat_1, x_hat_2, ...
%                                 f_nyquist, f_s1, f_s2, quantizer_name, M)
%     % PLOTQUANTIZERBYMVALUE Creates a single figure comparing all three frequencies for a specific M value
%     % Shows quantized signal and reconstruction side-by-side for each frequency
%     %
%     % Inputs:
%     %   t_n, x_n, xq_n - Nyquist sampled time, signal, and quantized signal
%     %   t_1, x_1, xq_1 - Undersampled time, signal, and quantized signal
%     %   t_2, x_2, xq_2 - Oversampled time, signal, and quantized signal
%     %   t_rec, xt_rec - Reconstruction time and original signal
%     %   x_hat_n, x_hat_1, x_hat_2 - Reconstructed signals for each frequency
%     %   f_nyquist, f_s1, f_s2 - Sampling frequencies
%     %   quantizer_name - Name of quantizer
%     %   M - Quantization levels
    
%     figure('Name', sprintf('%s M = %d', quantizer_name, M));
    
%     % Nyquist case: quantized (top-left) and reconstructed (top-right)
%     subplot(3, 2, 1);
%     plot(t_n, x_n, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_n, xq_n, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Nyquist (f_s = %.2f Hz): Quantized', f_nyquist), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 2, 2);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat_n, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Nyquist (f_s = %.2f Hz): Reconstructed', f_nyquist), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     % Undersampled case: quantized (middle-left) and reconstructed (middle-right)
%     subplot(3, 2, 3);
%     plot(t_1, x_1, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_1, xq_1, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Undersampled (f_s = %.2f Hz): Quantized', f_s1), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 2, 4);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat_1, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Undersampled (f_s = %.2f Hz): Reconstructed', f_s1), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     % Oversampled case: quantized (bottom-left) and reconstructed (bottom-right)
%     subplot(3, 2, 5);
%     plot(t_2, x_2, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_2, xq_2, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Oversampled (f_s = %.2f Hz): Quantized', f_s2), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 2, 6);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat_2, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Oversampled (f_s = %.2f Hz): Reconstructed', f_s2), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
% end

% function plotFrequencyComparison(t_samp, x_samp, xq_2lvl, xq_u1, xq_u2, xq_u3, xq_l1, xq_l2, xq_l3, ...
%                                   t_rec, xt_rec, x_hat, f_s, f_nyquist, freq_label, M1, M2, M3)
%     % PLOTFREQUENCYCOMPARISON Creates a single figure comparing all quantizer types for a specific frequency
%     % Shows quantized signal and reconstruction side-by-side for each quantizer type
%     %
%     % Inputs:
%     %   t_samp, x_samp - Sampled time and signal
%     %   xq_2lvl, xq_u1, xq_u2, xq_u3, xq_l1, xq_l2, xq_l3 - Quantized signals for each quantizer
%     %   t_rec, xt_rec, x_hat - Reconstruction time, original signal, and reconstructed signal
%     %   f_s - Sampling frequency
%     %   f_nyquist - Nyquist frequency
%     %   freq_label - Label for this frequency case
%     %   M1, M2, M3 - Quantization levels for uniform and Lloyd-Max
    
%     figure('Name', sprintf('Quantization: %s (f_s = %.2f Hz)', freq_label, f_s));
    
%     % Two-level quantizer (top row)
%     subplot(4, 2, 1);
%     plot(t_samp, x_samp, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_samp, xq_2lvl, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 3, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title('Two-level: Quantized', 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     subplot(4, 2, 2);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.2, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title('Two-level: Reconstructed', 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     % Uniform M1 (second row)
%     subplot(4, 2, 3);
%     plot(t_samp, x_samp, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_samp, xq_u1, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 3, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title(sprintf('Uniform M=%d: Quantized', M1), 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     subplot(4, 2, 4);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.2, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title(sprintf('Uniform M=%d: Reconstructed', M1), 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     % Uniform M2 (third row)
%     subplot(4, 2, 5);
%     plot(t_samp, x_samp, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_samp, xq_u2, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 3, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title(sprintf('Uniform M=%d: Quantized', M2), 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     subplot(4, 2, 6);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.2, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title(sprintf('Uniform M=%d: Reconstructed', M2), 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     % Uniform M3 and Lloyd-Max comparison (fourth row)
%     subplot(4, 2, 7);
%     plot(t_samp, x_samp, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_samp, xq_u3, 'DisplayName', 'Uniform', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 3, ...
%          'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     stem(t_samp, xq_l3, 'DisplayName', 'Lloyd-Max', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', 3, ...
%          'Color', [0.93 0.69 0.13], 'MarkerFaceColor', [0.93 0.69 0.13]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title(sprintf('M=%d Comparison: Quantized', M3), 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
    
%     subplot(4, 2, 8);
%     plot(t_rec, xt_rec, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec, x_hat, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.2, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('Amplitude', 'FontSize', 9);
%     title(sprintf('M=%d Comparison: Reconstructed', M3), 'FontSize', 10);
%     legend('show', 'Location', 'best', 'FontSize', 7);
%     grid on;
%     set(gca, 'FontSize', 8, 'LineWidth', 1);
% end

% function plotQuantizedAndReconstructed(t_n, x_n, xq_n, t_rec_n, xt_n, x_hat_n, ...
%                                         t_1, x_1, xq_1, t_rec_1, xt_1, x_hat_1, ...
%                                         t_2, x_2, xq_2, t_rec_2, xt_2, x_hat_2, ...
%                                         f_s1, f_s2, quantizer_name, M)
%     % PLOTQUANTIZEDANDRECONSTRUCTED Creates a single figure with quantized signal and reconstruction
%     % for three sampling cases, stacked vertically (2 subplots per case: quantized on left, reconstructed on right)
%     %
%     % Inputs:
%     %   t_n, x_n, xq_n - Nyquist sampled data and quantized version
%     %   t_rec_n, xt_n, x_hat_n - Reconstruction time, original signal, reconstructed signal (Nyquist)
%     %   t_1, x_1, xq_1 - Undersampled data and quantized version
%     %   t_rec_1, xt_1, x_hat_1 - Reconstruction for undersampled case
%     %   t_2, x_2, xq_2 - Oversampled data and quantized version
%     %   t_rec_2, xt_2, x_hat_2 - Reconstruction for oversampled case
%     %   f_s1, f_s2 - Sampling frequencies
%     %   quantizer_name - Name of quantizer
%     %   M - Number of quantization levels
    
%     figure('Name', sprintf('%s M = %d', quantizer_name, M));
    
%     % Nyquist case: quantized (top-left) and reconstructed (top-right)
%     subplot(3, 2, 1);
%     plot(t_n, x_n, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_n, xq_n, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title('Nyquist: Quantized', 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 2, 2);
%     plot(t_rec_n, xt_n, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec_n, x_hat_n, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title('Nyquist: Reconstructed', 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     % Undersampled case: quantized (middle-left) and reconstructed (middle-right)
%     subplot(3, 2, 3);
%     plot(t_1, x_1, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_1, xq_1, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Undersampled (f_s = %.2f Hz): Quantized', f_s1), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 2, 4);
%     plot(t_rec_1, xt_1, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec_1, x_hat_1, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Undersampled (f_s = %.2f Hz): Reconstructed', f_s1), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     % Oversampled case: quantized (bottom-left) and reconstructed (bottom-right)
%     subplot(3, 2, 5);
%     plot(t_2, x_2, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_2, xq_2, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Oversampled (f_s = %.2f Hz): Quantized', f_s2), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 2, 6);
%     plot(t_rec_2, xt_2, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t_rec_2, x_hat_2, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Oversampled (f_s = %.2f Hz): Reconstructed', f_s2), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
% end

% function plotQuantizedSignals(t_sample_nqyuist, x_sample_nyquist, xq_n, t_sample1, x_sample1, xq_1, t_sample2, x_sample2, xq_2, f_s1, f_s2, quantizer_name, M)
%     % PLOTQUANTIZEDSIGNALS Creates a 3-subplot figure showing quantized signals
%     %
%     % Inputs:
%     %   t_sample_nqyuist, x_sample_nyquist, xq_n - Nyquist sampled data
%     %   t_sample1, x_sample1, xq_1 - Aliased sampled data
%     %   t_sample2, x_sample2, xq_2 - Validly sampled data
%     %   f_s1, f_s2 - Sampling frequencies
%     %   quantizer_name - Name of quantizer (e.g., 'Uniform Quantization')
%     %   M - Number of quantization levels
    
%     figure('Name', sprintf('%s M = %d', quantizer_name, M));
    
%     subplot(3, 1, 1);
%     plot(t_sample_nqyuist, x_sample_nyquist, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_sample_nqyuist, xq_n, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title('Nyquist Sampling', 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 1, 2);
%     plot(t_sample1, x_sample1, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_sample1, xq_1, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 1, 3);
%     plot(t_sample2, x_sample2, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t_sample2, xq_2, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
% end

% function plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, quantizer_name, M)
%     % PLOTRECONSTRUCTEDSIGNALS Creates a 3-subplot figure showing reconstructed signals
%     %
%     % Inputs:
%     %   t, xt - Original continuous signal
%     %   x_hat_nyquist, x_hat_1, x_hat_2 - Reconstructed signals
%     %   f_s1, f_s2 - Sampling frequencies
%     %   quantizer_name - Name of quantizer
%     %   M - Number of quantization levels
    
%     figure('Name', sprintf('%s M = %d Reconstructed Signals', quantizer_name, M));
    
%     subplot(3, 1, 1);
%     plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t, x_hat_nyquist, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title('Nyquist Sampling', 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 1, 2);
%     plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t, x_hat_1, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
    
%     subplot(3, 1, 3);
%     plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     plot(t, x_hat_2, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 11);
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
% end

% function plotXq(t, x_samples, xq)
    % PLOTXQ Plots the original and quantized signals for visual comparison
    %
    % This function generates a plot that compares the original signal 
    % samples with the quantized signal samples produced by the Lloyd-Max 
    % quantization process. It helps in visualizing the effect of quantization 
    % on the signal.
    %
    % Inputs:
    %   t          - A vector of time samples corresponding to the signal.
    %   x_samples  - A vector of original signal samples.
    %   xq         - A vector of quantized signal samples.
    %
    % Example:
    %   plotXq(t, x_samples, xq);
%     plot(t, x_samples, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
%     hold on;
%     stem(t, xq, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
%     hold off;
%     legend('show', 'Location', 'best', 'FontSize', 8);
%     xlabel('Time (s)', 'FontSize', 10);
%     ylabel('Amplitude', 'FontSize', 10);
%     grid on;
%     set(gca, 'FontSize', 9, 'LineWidth', 1);
% end