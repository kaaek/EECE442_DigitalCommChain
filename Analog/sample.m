function [tSample, xSample] = sample(t, xt, F)
    % SAMPLE Samples continuous signal at specified sampling frequency
    %   Input: t (time vector), xt (signal values), F (sampling frequency Hz)
    %   Output: tSample (sample times), xSample (sampled values)
    T = 1/F;                                        % Define step size
    tSample = t(1):T:t(end);                        % Traverse from the signal's support in step size T    
    xSample = interp1(t, xt, tSample, 'spline');    % 'spline' option found by trial and error. Check interp1 documentation.
end