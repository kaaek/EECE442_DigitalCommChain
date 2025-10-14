% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [lvl, thr] = lloydMaxInit(x_samples, M)
    lvl = quantile(x_samples, (0.5:1:M)/M);     % Baseline case for representation points (equidistant)
    thr = (lvl(1:end-1)+lvl(2:end))/2;          % Define regions according to the selected representation points (with implicit t_0 = -inf and t_M = +inf)
end