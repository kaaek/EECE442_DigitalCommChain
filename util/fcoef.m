% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function cj = fcoef(t, xt, T, j) % Helper function that calculates the jth fourier series coefficient.
    %FCOEF   Fourier series coefficient calculation
    %   CJ = FCOEF(T, XT, T, J) computes the j-th Fourier series coefficient
    %   for the signal XT defined on the time axis T with period T.
    %
    %   Inputs:
    %       T  - row vector, time axis
    %       XT - row vector, signal values
    %       T  - scalar, period of the signal
    %       J  - scalar, index of the coefficient to compute
    %
    %   Output:
    %       CJ - scalar, j-th Fourier series coefficient
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       T = 1/10; % Period of the signal
    %       j = 1;
    %       cj = fcoef(t, xt, T, j);
    %       disp(cj);
    cj = (1/T) * trapz(t, xt .* exp(-1i*2*pi*j*t/T));
end