function xhat = fs(xt, t, N, T)
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
    
    % Initialize
    k = -N:N; % Total number of coefficients to compute is 2*N + 1 (symmetry and zero)
    ck = zeros(1, length(k));
    xhat = zeros(size(t));
    for i = 1:length(k)
        ck(i) = fcoef(t, xt, T, k(i));                  % Calculate coefficients c_k
        xhat = xhat + ck(i) * exp(1i*2*pi*k(i)*t/T);    % Reconstruct signal using c_k
    end
    xhat = real(xhat); % So that we can analyze practically & plot.
end