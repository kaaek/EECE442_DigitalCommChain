% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

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
    %           By definition, this is: t = n*T_s
    %       X_SAMPLE - row vector, sampled signal values.
    %           By definition this is: x[n] = x(t) at t = nT_s such that n is in the integer set Z.
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       [ts, xs] = sample(t, xt, 20);
    %       plot(ts, xs, 'o')
    Ts = 1/fs;
    t_sample = t(1):Ts:t(end);
    x_sample = interp1(t, xt, t_sample, 'spline'); % 'spline' option found by trial and error. Check interp1 documentation.

    figure;
    plot(t, xt, 'r--', 'DisplayName', 'x(t)'); hold on;
    stem(t_sample, x_sample, 'b', 'DisplayName', 'x[n]', 'Marker', 'o');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    title(sprintf('[TX] Sampler: x(t) vs x[n] sampled at f_s = %.2f Hz', fs));
    legend show;
    grid on;
end