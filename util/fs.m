% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function xhat = fs(xt, t, n, T)
    % FS   Fourier Series approximation
    %   XHAT = FS(XT, T, N, T) computes the Fourier Series
    %   approximation of the signal XT defined on the time axis T using
    %   N coefficients.
    %
    %   Inputs:
    %       XT - row vector, signal values
    %       T  - row vector, time axis
    %       N  - scalar, number of coefficients to use
    %       T  - scalar, period of the signal
    %
    %   Outputs:
    %       XHAT - row vector, approximated signal values
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       n = 5;
    %       T = 1/10; % Period of the signal
    %       [xhat, ck] = fs(xt, t, n, T);
    %       plot(t, xt, 'r-', t, xhat, 'b--')
    %       legend('Original Signal', 'Fourier Series Approximation');
    k = -n:n;
    ck = zeros(1, length(k));
    xhat = zeros(size(t));

    for idx = 1:length(k)   % Calculate coefficients c_k
        ck(idx) = fcoef(t, xt, T, k(idx));
    end
    for idx = 1:length(k)   % Reconstruct signal using c_k
        xhat = xhat + ck(idx) * exp(1i*2*pi*k(idx)*t/T);
    end

    
    plot(t, xt, xhat);


end

function plot(t, xt, xhat)
    figure;
    plot(t, xt, 'r--', 'LineWidth', 1.5);
    hold on;
    plot(t, real(xhat), 'b-', 'LineWidth', 1.5);
    hold off;
    legend('x(t)', 'x^(t)');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    title('[TX] Fourier Series: x(t) & x^(t) VS Time (s)');
    grid on;
end