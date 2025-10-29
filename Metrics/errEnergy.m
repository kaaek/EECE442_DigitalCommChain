% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function E = errEnergy (xt, xhat, t)
    % ERREREN   Error Energy Calculation
    %   E = ERREREN(T, XT, XHAT) computes the error energy between the 
    %   original signal XT and the approximated signal XHAT over the time 
    %   axis T.
    %
    %   Inputs:
    %       T    - row vector, time axis
    %       XT   - row vector, original signal values
    %       XHAT - row vector, approximated signal values
    %
    %   Output:
    %       E - scalar, computed error energy
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       xhat = xt + 0.1*randn(size(xt)); % Example of an approximated signal
    %       E = errEn(t, xt, xhat);
    
    integrand = (xt - xhat).^2;
    E = integrate(t, integrand);
end