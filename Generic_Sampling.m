clc, clearvars % Flush everything before running

function [t_sample, x_sample] = sample(t, xt, fs)
    %SAMPLE   Sample a signal at a given frequency
    %   [T_SAMPLE, X_SAMPLE] = SAMPLE(T, XT, FS) resamples the input signal
    %   XT defined on the time axis T at a sampling frequency FS.
    %
    %   Inputs:
    %       T   - row vector, time axis
    %       XT  - row vector, signal values
    %       FS  - scalar, sampling frequency
    %
    %   Outputs:
    %       T_SAMPLE - row vector, sampled time axis
    %       X_SAMPLE - row vector, sampled signal values
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       [ts, xs] = sample(t, xt, 20);
    %       plot(ts, xs, 'o')
    t_sample = t(1):1/fs:t(end);
    x_sample = interp1(t, xt, t_sample, 'spline'); % 'spline' option found by trial and error. Check interp1 documentation.
end

function y = sinc(x)
    y = ones(size(x));
    idx = x ~= 0;
    y(idx) = sin(pi*x(idx))./(pi*x(idx));
end


function xrcon = reconstruct(t, x_sample, fs)
    % Help %
    Ts = 1/fs;
    n = 0:length(x_sample)-1;     % sample indices
    t_sample = n*Ts;              % sample times
    xrcon = zeros(size(t));
    for k = 1:length(n)
        xrcon = xrcon + x_sample(k) * sinc(fs*(t - t_sample(k)));
    end
end

function [xhat, ck] = ffs(xt, t, n, T)
    k = -n:n;
    ck = zeros(1, length(k));
    xhat = zeros(size(t));

    % Compute coefficients
    for idx = 1:length(k)
        ck(idx) = fcoef(t, xt, T, k(idx));
    end

    % Reconstruct signal
    for idx = 1:length(k)
        xhat = xhat + ck(idx) * exp(1i*2*pi*k(idx)*t/T);
    end
end

function cj = fcoef(t, xt, T, j)
    cj = (1/T) * trapz(t, xt .* exp(-1i*2*pi*j*t/T));
end