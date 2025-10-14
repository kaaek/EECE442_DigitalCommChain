function [n, xrcon] = reconstruct(t, x_sample, fs)
    %RECONSTRUCT   Reconstruct a signal from its samples
    %   XRCON = RECONSTRUCT(T, X_SAMPLE, FS) reconstructs the original signal
    %   from the sampled signal X_SAMPLE defined on the time axis T at a 
    %   sampling frequency FS.
    %
    %   Inputs:
    %       T        - row vector, time axis for the original signal
    %       X_SAMPLE - row vector, sampled signal values
    %       FS       - scalar, sampling frequency used for sampling
    %
    %   Output:
    %       XRCON    - row vector, reconstructed signal values
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       fs = 20;
    %       [ts, xs] = sample(t, xt, fs);
    %       xr = reconstruct(t, xs, fs);
    %       plot(t, xt, 'r-', t, xr, 'b--')
    %       legend('Original Signal', 'Reconstructed Signal');
    Ts = 1/fs;
    n = 0:length(x_sample)-1;
    t_sample = n*Ts;              
    xrcon = zeros(size(t));
    for k = 1:length(n)
        xrcon = xrcon + x_sample(k) * sinc(fs * (t - t_sample(k)));
    end
end