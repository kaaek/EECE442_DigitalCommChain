% ----------------------------------------------------------------------
% authors: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function MSE = MSE (t, xt, xhat)
    % MSE calculates the Mean Squared Error between the true signal and the estimated signal.
    % 
    % Inputs:
    %   t      - Time vector (1D array)
    %   xt     - True signal (1D array)
    %   xhat   - Estimated signal (1D array)
    %
    % Output:
    %   MSE    - Mean Squared Error (scalar)
    %
    % Usage:
    %   mseValue = MSE(t, xt, xhat);
    MSE = errEn(t, xt, xhat)/length(t);
end