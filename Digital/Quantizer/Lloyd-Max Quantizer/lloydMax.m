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
    xq = quantize(x_samples, thr, lvl);
    MSE = mean((x_samples - xq).^2);
    % fprintf('\nLLOYD-MAX QUANTIZER (M = %s), ', num2str(M));
    % fprintf('MEAN SQUARE ERROR (MSE): %.6f\n', MSE);
end

function f_x = map(x, thr, lvl)
    % MAP Maps input samples to quantized levels based on thresholds.
    %
    % This function takes a vector of input samples and maps them to 
    % the corresponding quantization levels based on the provided 
    % thresholds. The mapping is done by checking which threshold 
    % interval each sample falls into.
    %
    % Inputs:
    %   x   - A vector of input samples to be quantized.
    %   thr - A vector of thresholds used for quantization.
    %   lvl - A vector of representation levels corresponding to the thresholds.
    %
    % Outputs:
    %   f_x - A vector of quantized output samples after mapping.
    
    f_x = zeros(size(x)); % Initialize the output vector with zeros.
    for i = 1:length(x)
        if x(i) < thr(1) % If the sample is less than the first threshold
            f_x(i) = lvl(1); % Assign the first level
        elseif x(i) >= thr(end) % If the sample is greater than or equal to the last threshold
            f_x(i) = lvl(end); % Assign the last level
        else
            for j = 1:length(thr)-1 % Iterate through the thresholds
                if x(i) >= thr(j) && x(i) < thr(j+1) % Check which interval the sample falls into
                    f_x(i) = lvl(j+1); % Assign the corresponding level
                    break; % Exit the loop once the level is assigned
                end
            end
        end
    end
end

function mse = calcMSE(x, f_x)
    mse = mean(abs(x - f_x).^2); % By the law of large numbers
end