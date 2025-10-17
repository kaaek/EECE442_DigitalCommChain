% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [xq, MSE] = uniformQuan(M, t, xt, print)
    % uniformQuan quantizes a signal using uniform quantization.
    %
    % Syntax:
    %   [xq, MSE] = uniformQuan(M, t, xt)
    %
    % Inputs:
    %   M   - A vector containing the number of quantization levels.
    %   t   - A vector representing the time instances corresponding to the signal.
    %   xt  - The original signal to be quantized.
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
    % MSE = zeros(1, length(M));
    % figure('Name','Uniform Quantizer');
    % for i = 1:length(M) % Plot the quantized signal
    lvl = linspace(min(xt), max(xt), M);
    thr = (lvl(1:end-1) + lvl(2:end)) / 2;
    xq = quan(xt, thr, lvl);
    MSE = mean((xt - xq).^2);
    
    % Print the quantization levels and thresholds
    
    if print
        fprintf('### Number of Quantization Levels: %s ###\n', num2str(M));
        fprintf('Quantization Levels:\n');
        fprintf('%.4f ', lvl);
        fprintf('\nThresholds:\n');
        fprintf('%.4f ', thr);
        fprintf('\n\n\n');
    end
end