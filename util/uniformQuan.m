% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [thr, lvl, xq, MSE] = uniformQuan(M, t, xt, print)
    % uniformQuan quantizes a signal using uniform quantization.
    %
    % Syntax:
    %   [xq, MSE] = uniformQuan(M, t, xt)
    %
    % Inputs:
    %   M   - A scalar containing the number of quantization levels.
    %   t   - A vector representing the time instances corresponding to the signal.
    %   xt  - The original signal to be quantized (vector).
    %
    % Outputs:
    %   xq  - The quantized signal.
    %   MSE - A vector containing the Mean Squared Error for each quantization level.
    %
    % Description:
    %   This function performs uniform quantization on the input signal 'xt'
    %   for each value in the vector 'M'. It computes the quantized signal
    %   and the Mean Squared Error (MSE) between the original and quantized
    %   signals. The function also generates plots to visualize the original
    %   and quantized signals, as well as the MSE against the quantization levels.
    lvl = linspace(min(xt), max(xt), M);
    thr = (lvl(1:end-1) + lvl(2:end)) / 2; % Calculate the quantization levels and thresholds for uniform quantization
    xq = quan(xt, thr, lvl); % Quantize the original signal using the calculated thresholds and levels
    MSE = mean((xt - xq).^2);
    
    if print
        fprintf('-----------------------------------------------------------------');
        fprintf('\nUNIFORM QUANTIZER (Number of Quantization Levels M = %s)', num2str(M));
        fprintf('\nQUANTIZATION LEVELS:\n');
        fprintf('%.4f ', lvl);
        fprintf('\nTHRESHOLDS:\n');
        fprintf('%.4f ', thr);
        fprintf('\nMEAN SQUARED ERROR (MSE): %.4e\n', MSE);
        fprintf('\n\n\n');
    end
end