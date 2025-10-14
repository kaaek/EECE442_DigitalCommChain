function xhat = fs(xt, t, n, T)
    % FS   Fourier Series approximation
    %   [XHAT, CK] = FS(XT, T, N, T) computes the Fourier Series
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
    %       CK   - row vector, Fourier coefficients
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
end