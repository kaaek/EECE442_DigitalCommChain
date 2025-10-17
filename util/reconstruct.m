% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function xrcon = reconstruct(t_target, x_sample, fs)
    % RECONSTRUCT reconstructs a continuous signal from its sampled values.
    %
    %   XRCON = RECONSTRUCT(T_TARGET, X_SAMPLE, FS) takes a target time vector
    %   T_TARGET, a vector of sampled signal values X_SAMPLE, and the sampling
    %   frequency FS. It reconstructs the continuous signal XRCON at the time
    %   points specified in T_TARGET using the sinc interpolation method.
    %
    %   Inputs:
    %       T_TARGET - A vector of time points at which to reconstruct the signal.
    %       X_SAMPLE - A vector of sampled signal values.
    %       FS       - The sampling frequency of the original signal.
    %
    %   Output:
    %       XRCON    - A vector containing the reconstructed signal values at
    %                   the specified time points in T_TARGET.
    t_sample = (0:length(x_sample)-1)/fs;
    xrcon = zeros(size(t_target));
    for k = 1:length(x_sample)
        xrcon = xrcon + x_sample(k) * sinc(fs*(t_target - t_sample(k)));
    end
end
