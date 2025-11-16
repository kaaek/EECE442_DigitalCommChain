% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function E = errorEnergy (xt, xhat, t)
    integrand = (xt - xhat).^2;
    E = integrate(t, integrand);
end