% ----------------------------------------------------------------------
% authors: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function foo = MSE (xt, xhat, t)
    foo = errorEnergy(t, xt, xhat)/length(t);
end