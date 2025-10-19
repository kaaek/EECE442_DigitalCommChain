% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function y = sinc(x)
    % SINC Compute the normalized sinc function.
    %
    %   Y = SINC(X) computes the normalized sinc function defined as:
    %   sinc(x) = sin(pi*x) / (pi*x) for x ~= 0, and sinc(0) = 1.
    %
    %   Input:
    %       X - A scalar, vector, or matrix of values at which to evaluate the sinc function.
    %
    %   Output:
    %       Y - The computed sinc values, returned in the same shape as X.
    y = ones(size(x));
    idx = x ~= 0;
    y(idx) = sin(pi*x(idx))./(pi*x(idx));
end