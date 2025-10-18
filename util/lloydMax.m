% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [thr, lvl, xq, MSE] = lloydMax(x_samples, M, tgtMSE)
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
    
    [lvl, thr] = lloydMaxInit(x_samples, M);                % Initialize with equidistant representation points.
    iter = 0;
    maxIter = inf;                                          % Used as a failsafe limit, can make this finite if needed.
    currMSE = inf;
    while currMSE > tgtMSE && iter < maxIter
        iter = iter + 1;
        lvl_prev = lvl;
        thr = (lvl(1:end-1) + lvl(2:end)) / 2;              % Update thresholds based on current levels (mean)
        lvl = partition(x_samples, thr, M, lvl_prev);       % Partition the samples and update levels
        f_x = map(x_samples, thr, lvl);
        currMSE = calcMSE(x_samples, f_x);                  % Calculate the current MSE. Check if it falls below the target.
    end
    xq = quan(x_samples, thr, lvl);
    MSE = mean((x_samples - xq).^2);
    fprintf('-----------------------------------------------------------------');
    fprintf('\nLLOYD-MAX QUANTIZER (Number of Quantization Levels M = %s)', num2str(M));
    fprintf('\nMEAN SQUARE ERROR (MSE): %.4f\n', MSE);
    fprintf('\n\n\n');
end

function f_x = map(x, thr, lvl)
    f_x = zeros(size(x));
    for i = 1:length(x)
        if x(i) < thr(1)
            f_x(i) = lvl(1);
        elseif x(i) >= thr(end)
            f_x(i) = lvl(end);
        else
            for j = 1:length(thr)-1
                if x(i) >= thr(j) && x(i) < thr(j+1)
                    f_x(i) = lvl(j+1);
                    break;
                end
            end
        end
    end
end

function mse = calcMSE(x, f_x)
    mse = mean(abs(x - f_x).^2); % By the law of large numbers
end