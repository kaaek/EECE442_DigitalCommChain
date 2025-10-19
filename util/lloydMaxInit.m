% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [lvl, thr] = lloydMaxInit(x_samples, M)
    % LLOYDMAXINIT initializes the representation levels and thresholds for the 
    % Lloyd-Max quantization algorithm. This function serves as a helper method 
    % for the main Lloyd-Max function (lloydMax.m). It computes the initial 
    % representation points (levels) based on the provided sample data and 
    % the number of quantization levels (M). The thresholds are then calculated 
    % to define the regions for quantization.
    % 
    % Inputs:
    %   x_samples - A vector of sample data points to be quantized.
    %   M         - The number of quantization levels.
    %
    % Outputs:
    %   lvl - A vector containing the initial representation levels.
    %   thr - A vector containing the thresholds that define the quantization 
    %         regions.
    % lvl = quantile(x_samples, (0.5:1:M)/M);     % Baseline case for representation points (equidistant)
    lvl = linspace(min(x_samples), max(x_samples), M);
    thr = (lvl(1:end-1)+lvl(2:end))/2;          % Define regions according to the selected representation points (with implicit t_0 = -inf and t_M = +inf)
end