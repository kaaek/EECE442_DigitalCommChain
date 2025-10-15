% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function xq = quan(x, thr, lvl)
    % quan - Quantizes the input signal based on specified thresholds and levels.
    %
    % Syntax:
    %   xq = quan(x, thr, lvl)
    %
    % Inputs:
    %   x   - Input signal to be quantized mse vs quantization level M(vector).
    %   thr - Thresholds for quantization (vector). The length of thr should be one less than the length of lvl.
    %   lvl - Levels corresponding to the thresholds (vector). The length of lvl should be one more than the length of thr.
    %
    % Outputs:
    %   xq  - Quantized output signal (vector) of the same size as input x.
    %
    % Description:
    %   The function quantizes the input signal x based on the provided thresholds and levels. 
    %   Each value in x is replaced by the corresponding level based on the defined thresholds.
    %   If a value in x is less than or equal to the first threshold, it is assigned the first level.
    %   If it is greater than the last threshold, it is assigned the last level. 
    %   Values between thresholds are assigned the corresponding levels.
    xq = zeros(size(x));
    for i = 1:length(x)
        x_i = x(i);
        if x_i <= thr(1)
            xq(i) = lvl(1);
        elseif x_i > thr(end)
            xq(i) = lvl(end);
        else
            for j = 1:length(thr)-1
                if x_i > thr(j) && x_i <= thr(j+1)
                    xq(i) = lvl(j+1);
                    break;
                end
            end
        end
    end
end