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

% The two-level quantizer
xqn = twoLvlQuan(x_sample_nyquist);
xq1 = twoLvlQuan(x_sample1);
xq2 = twoLvlQuan(x_sample2);
plotQuantizedSignals(t_sample_nqyuist, x_sample_nyquist, xqn, t_sample1, x_sample1, xq1, t_sample2, x_sample2, xq2, f_s1, f_s2, 'Two Level Quantizer', 2);
% Reconstruct an approximate signal x^(t) from the quantized signal via sinc-interpolation
x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2  = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Two-level Quantizer', 2);

% ----------------------------------------------------------------------

M1 = 4;
M2 = 16;
M3 = 32;

[xqn, ~, ~, ~] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M1);
[xq1, ~, ~, ~] = uniformQuan(t_sample1, x_sample1, M1);
[xq2, ~, ~, ~] = uniformQuan(t_sample2, x_sample2, M1);

plotQuantizedSignals(t_sample_nqyuist, x_sample_nyquist, xqn, t_sample1, x_sample1, xq1, t_sample2, x_sample2, xq2, f_s1, f_s2, 'Uniform Quantization', M1);

x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2 = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Uniform Quantizer', M1);

% ----------------------------------------------------------------------

[xqn, ~, ~, ~] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M2);
[xq1, ~, ~, ~] = uniformQuan(t_sample1, x_sample1, M2);
[xq2, ~, ~, ~] = uniformQuan(t_sample2, x_sample2, M2);

plotQuantizedSignals(t_sample_nqyuist, x_sample_nyquist, xqn, t_sample1, x_sample1, xq1, t_sample2, x_sample2, xq2, f_s1, f_s2, 'Uniform Quantization', M2);

x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2 = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Uniform Quantizer', M2);

% -------------------------------------------------------------------------

[xqn, ~, ~, ~] = uniformQuan(t_sample_nqyuist, x_sample_nyquist, M3);
[xq1, ~, ~, ~] = uniformQuan(t_sample1, x_sample1, M3);
[xq2, ~, ~, ~] = uniformQuan(t_sample2, x_sample2, M3);

plotQuantizedSignals(t_sample_nqyuist, x_sample_nyquist, xqn, t_sample1, x_sample1, xq1, t_sample2, x_sample2, xq2, f_s1, f_s2, 'Uniform Quantization', M3);

x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2 = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Uniform Quantizer', M3);

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
xlabel('Number of Quantization Levels (M)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean Squared Error (MSE)', 'FontSize', 12, 'FontWeight', 'bold');
title('MSE vs Number of Quantization Levels', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);

% Near-optimal quantization and overall assessment
tgtMSE = 0.1;

[thr_n, lvl_n, xq_n, ~] = lloydMax(x_sample_nyquist, M1, tgtMSE);
[thr_1, lvl_1, xq_1, ~] = lloydMax(x_sample1, M1, tgtMSE);
[thr_2, lvl_2, xq_2, ~] = lloydMax(x_sample2, M1, tgtMSE);

figure('Name',sprintf('Lloyd-Max Quantized Signals M = %d', M1));
subplot(3, 1, 1);
plotXq(t_sample_nqyuist, x_sample_nyquist, xq_n);
title('Nyquist Sampling', 'FontSize', 12, 'FontWeight', 'bold');

subplot(3, 1, 2);
plotXq(t_sample1, x_sample1, xq_1);
title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 12, 'FontWeight', 'bold');

subplot(3, 1, 3);
plotXq(t_sample2, x_sample2, xq_2);
title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 12, 'FontWeight', 'bold');

plotPdf(x_sample_nyquist, thr_n, lvl_n, M1);
plotPdf(x_sample1, thr_1, lvl_1, M1);
plotPdf(x_sample2, thr_2, lvl_2, M1);

x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2 = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Lloyd-Max Quantizer', M1);

% --------------------------------------------------------------------------------------------------

[thr_n, lvl_n, xq_n, ~] = lloydMax(x_sample_nyquist, M2, tgtMSE);
[thr_1, lvl_1, xq_1, ~] = lloydMax(x_sample1, M2, tgtMSE);
[thr_2, lvl_2, xq_2, ~] = lloydMax(x_sample2, M2, tgtMSE);

figure('Name',sprintf('Lloyd-Max Quantized Signals M = %d', M2));
subplot(3, 1, 1);
plotXq(t_sample_nqyuist, x_sample_nyquist, xq_n);
title('Nyquist Sampling', 'FontSize', 12, 'FontWeight', 'bold');

subplot(3, 1, 2);
plotXq(t_sample1, x_sample1, xq_1);
title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 12, 'FontWeight', 'bold');

subplot(3, 1, 3);
plotXq(t_sample2, x_sample2, xq_2);
title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 12, 'FontWeight', 'bold');

plotPdf(x_sample_nyquist, thr_n, lvl_n, M2);
plotPdf(x_sample1, thr_1, lvl_1, M2);
plotPdf(x_sample2, thr_2, lvl_2, M2);

x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2 = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Lloyd-Max Quantizer', M2);

% --------------------------------------------------------------------------------------------------

[thr_n, lvl_n, xq_n, ~] = lloydMax(x_sample_nyquist, M3, tgtMSE);
[thr_1, lvl_1, xq_1, ~] = lloydMax(x_sample1, M3, tgtMSE);
[thr_2, lvl_2, xq_2, ~] = lloydMax(x_sample2, M3, tgtMSE);

figure('Name',sprintf('Lloyd-Max Quantized Signals M = %d', M3));
subplot(3, 1, 1);
plotXq(t_sample_nqyuist, x_sample_nyquist, xq_n);
title('Nyquist Sampling', 'FontSize', 12, 'FontWeight', 'bold');

subplot(3, 1, 2);
plotXq(t_sample1, x_sample1, xq_1);
title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 12, 'FontWeight', 'bold');

subplot(3, 1, 3);
plotXq(t_sample2, x_sample2, xq_2);
title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 12, 'FontWeight', 'bold');

plotPdf(x_sample_nyquist, thr_n, lvl_n, M3);
plotPdf(x_sample1, thr_1, lvl_1, M3);
plotPdf(x_sample2, thr_2, lvl_2, M3);

x_hat_nyquist = reconstruct(t, x_sample_nyquist, f_Nyquist);
x_hat_1 = reconstruct(t, x_sample1, f_s1);
x_hat_2  = reconstruct(t, x_sample2, f_s2);

plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, 'Lloyd-Max Quantizer', M3);


end

function plotQuantizedSignals(t_sample_nqyuist, x_sample_nyquist, xq_n, t_sample1, x_sample1, xq_1, t_sample2, x_sample2, xq_2, f_s1, f_s2, quantizer_name, M)
    % PLOTQUANTIZEDSIGNALS Creates a 3-subplot figure showing quantized signals
    %
    % Inputs:
    %   t_sample_nqyuist, x_sample_nyquist, xq_n - Nyquist sampled data
    %   t_sample1, x_sample1, xq_1 - Aliased sampled data
    %   t_sample2, x_sample2, xq_2 - Validly sampled data
    %   f_s1, f_s2 - Sampling frequencies
    %   quantizer_name - Name of quantizer (e.g., 'Uniform Quantization')
    %   M - Number of quantization levels
    
    figure('Name', sprintf('%s M = %d', quantizer_name, M));
    
    subplot(3, 1, 1);
    plot(t_sample_nqyuist, x_sample_nyquist, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
    hold on;
    stem(t_sample_nqyuist, xq_n, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Nyquist Sampling', 'FontSize', 12, 'FontWeight', 'bold');
    legend('show', 'Location', 'best', 'FontSize', 9);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
    
    subplot(3, 1, 2);
    plot(t_sample1, x_sample1, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
    hold on;
    stem(t_sample1, xq_1, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 12, 'FontWeight', 'bold');
    legend('show', 'Location', 'best', 'FontSize', 9);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
    
    subplot(3, 1, 3);
    plot(t_sample2, x_sample2, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
    hold on;
    stem(t_sample2, xq_2, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 12, 'FontWeight', 'bold');
    legend('show', 'Location', 'best', 'FontSize', 9);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
end

function plotReconstructedSignals(t, xt, x_hat_nyquist, x_hat_1, x_hat_2, f_s1, f_s2, quantizer_name, M)
    % PLOTRECONSTRUCTEDSIGNALS Creates a 3-subplot figure showing reconstructed signals
    %
    % Inputs:
    %   t, xt - Original continuous signal
    %   x_hat_nyquist, x_hat_1, x_hat_2 - Reconstructed signals
    %   f_s1, f_s2 - Sampling frequencies
    %   quantizer_name - Name of quantizer
    %   M - Number of quantization levels
    
    figure('Name', sprintf('%s M = %d Reconstructed Signals', quantizer_name, M));
    
    subplot(3, 1, 1);
    plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
    hold on;
    plot(t, x_hat_nyquist, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Nyquist Sampling', 'FontSize', 12, 'FontWeight', 'bold');
    legend('show', 'Location', 'best', 'FontSize', 9);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
    
    subplot(3, 1, 2);
    plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
    hold on;
    plot(t, x_hat_1, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Undersampled (f_s = %.2f Hz)', f_s1), 'FontSize', 12, 'FontWeight', 'bold');
    legend('show', 'Location', 'best', 'FontSize', 9);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
    
    subplot(3, 1, 3);
    plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
    hold on;
    plot(t, x_hat_2, '-', 'DisplayName', 'Reconstructed', 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Oversampled (f_s = %.2f Hz)', f_s2), 'FontSize', 12, 'FontWeight', 'bold');
    legend('show', 'Location', 'best', 'FontSize', 9);
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
end

function plotXq(t, x_samples, xq)
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
    plot(t, x_samples, '--', 'DisplayName', 'Original', 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
    hold on;
    stem(t, xq, 'DisplayName', 'Quantized', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    hold off;
    legend('show', 'Location', 'best', 'FontSize', 9);
    xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 9, 'LineWidth', 1);
end

function plotPdf(x_samples, thr, lvl, M)
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
    %
    % Example:
    %   plotPdf(x_samples, thr, lvl);
    
    % Get the empirical PDF from the sampled sequence x_samples
    x_samples = real(x_samples);  % Disregard imaginary part
    nbins = min(max(80, round(length(x_samples) / 10)), 200);
    [pdf, edges] = histcounts(x_samples, nbins, 'Normalization', 'pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    figure('Name', ['Lloyd-Max Empirical PDF, M = ' num2str(M)]);
    plot(centers, pdf, 'LineWidth', 2, 'Color', [0 0.45 0.74]);
    hold on;
    for k = 1:length(thr)
        plot([thr(k) thr(k)], ylim, '--', 'LineWidth', 2, 'Color', [0.93 0.69 0.13]);
    end
    stem(lvl, interp1(centers, pdf, lvl, 'linear', 'extrap'), 'filled', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1]);
    legend('PDF', 'Thresholds', 'Levels', 'Location', 'best', 'FontSize', 10);
    title(sprintf('Lloyd-Max Empirical PDF (M = %d)', M), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end