% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [lvl, thr] = lloydMaxInit(x_samples, M)
    lvl = linspace(min(x_samples), max(x_samples), M);
    thr = (lvl(1:end-1)+lvl(2:end))/2;          % Define regions according to the selected representation points (with implicit t_0 = -inf and t_M = +inf)
end