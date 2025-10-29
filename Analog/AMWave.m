% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [t, xt, f_max] = AMWave(duration, FM, FC)
    t = 0:0.0001:duration;
    FM = 5;    % Message frequency
    f_max = FC + FM;
    xt = (1 + 0.5*cos(2*pi*FM*t)) .* cos(2*pi*FC*t); % Just an AM-modulated wave
end