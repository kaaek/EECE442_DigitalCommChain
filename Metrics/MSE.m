% ----------------------------------------------------------------------
% authors: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function MSE = MSE (xt, xhat, t)
    MSE = errEnergy(t, xt, xhat)/length(t);
end