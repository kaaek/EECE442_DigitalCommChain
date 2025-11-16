% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors if exists.
% ----------------------------------------------------------------------

function [tSample, xSample] = sample(t, xt, F)
    T = 1/F;                                        % Define step size
    tSample = t(1):T:t(end);                        % Traverse from the signal's support in step size T    
    xSample = interp1(t, xt, tSample, 'spline');    % 'spline' option found by trial and error. Check interp1 documentation.
end