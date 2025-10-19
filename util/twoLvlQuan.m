% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function xq = twoLvlQuan(t, xt)
    % twoLvlQuan - Performs two-level quantization on the input signal.
    %
    % Syntax: xq = twoLvlQuan(t, xt)
    %
    % Inputs:
    %   t  - Threshold value for quantization.
    %   xt - Input signal vector to be quantized.
    %
    % Outputs:
    %   xq - Quantized output signal vector.
    %
    % Description:
    %   This function takes an input signal and applies a two-level quantization
    %   based on the provided threshold. It calculates the mean of the input signal
    %   to determine the quantization levels and returns the quantized signal.
    t_1 = mean(xt);
    l_1 = (min(xt)+t_1)/2;
    l_2 = (max(xt)+t_1)/2;
    
    thr = [t_1];
    lvl = [l_1, l_2];
    xq = quan(xt, thr, lvl);
end