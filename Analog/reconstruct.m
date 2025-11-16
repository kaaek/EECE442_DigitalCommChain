% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function xrcon = reconstruct(T_TARGET, X_SAMPLE, F_S)
    t_sample = (0:length(X_SAMPLE)-1)/F_S;
    xrcon = zeros(size(T_TARGET));
    for k = 1:length(X_SAMPLE)
        xrcon = xrcon + X_SAMPLE(k) * sinc(F_S*(T_TARGET - t_sample(k))); % Add the contribution of the k-th sample to the reconstructed signal (with sinc interpolation)
    end
end
