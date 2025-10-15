% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [t, xt, f_max] = exampleSpeechWave(duration, f_c)
    t = 0:0.0001:duration;
    f_m = 5;    % Message frequency
    f_max = f_c + f_m;
    xt = (1 + 0.5*cos(2*pi*f_m*t)) .* cos(2*pi*f_c*t);
end