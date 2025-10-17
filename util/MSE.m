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