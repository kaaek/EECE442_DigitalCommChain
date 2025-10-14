% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [thr, lvl, xq, MSE] = lloydMax(x_samples, M)
    % LLOYDMAX Quantization using Lloyd-Max Algorithm
    %
    % This function implements the Lloyd-Max quantization algorithm to 
    % quantize a set of samples into M levels. The algorithm iteratively 
    % refines the representation levels and thresholds to minimize the 
    % quantization error.
    %
    % Inputs:
    %   x_samples - A vector of input samples to be quantized.
    %   M         - The number of quantization levels.
    %
    % Outputs:
    %   thr - A vector of thresholds used for quantization.
    %   lvl - A vector of representation levels after quantization.
    %   xq  - A vector of quantized output samples.
    %
    % Example:
    %   [thr, lvl, xq] = lloydMax(x_samples, M);
    
    [lvl, thr] = lloydMaxInit(x_samples, M);                % Initialize with equidistant representation points
    tol = 1e-8;                                             % Tolerance "epsilon"
    maxIter = 300;                                          % Maximum number of iterations (fail-safe)
    iter = 0;
    max_change = inf;                                       % Initialize to a large value

    while max_change >= tol && iter < maxIter               % Repeat until max_change < tol or iter >= maxIter
        iter = iter + 1;
        lvl_prev = lvl;
        thr = (lvl(1:end-1) + lvl(2:end)) / 2;              % Update thresholds based on current levels (mean)
        lvl = partition(x_samples, thr, M, lvl_prev);       % Partition the samples and update levels
        max_change = max(abs(lvl - lvl_prev));              % Calculate the maximum change in levels
    end
    xq = quan(x_samples, thr, lvl);
    MSE = mean((x_samples - xq).^2);
    plotPdf(x_samples, thr, lvl);
end

function plotPdf(x_samples, thr, lvl)
                                                                            % Get the empirical PDF from the sampled sequence x_samples
    nbins = min(max(80, round(length(x_samples) / 10)), 200);               % Determine number of bins based on data length
    [pdf, edges] = histcounts(x_samples, nbins, 'Normalization', 'pdf');    % Counts the number of occurences of the values in the sample → empirical PDF
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    figure;
    plot(centers, pdf, 'LineWidth', 2);
    hold on;
    for k = 1:length(thr)
        plot([thr(k) thr(k)], ylim, '--w');
    end
    stem(lvl, interp1(centers, pdf, lvl, 'linear', 'extrap'), 'r', 'filled');
    legend('Empirical PDF', 'Thresholds', 'Quantization Levels (red stems)');
    title('Sample Points PDF Curve, Representation Points, and Thresholds');
    xlabel('x');
    ylabel('Pr\{x\}');
end