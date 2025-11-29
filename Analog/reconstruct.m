function xrcon = reconstruct(T_TARGET, X_SAMPLE, F_S)
    % RECONSTRUCT Reconstructs continuous signal from samples using sinc interpolation
    %   Input: T_TARGET (target time points), X_SAMPLE (sample values), F_S (sampling rate Hz)
    %   Output: xrcon (reconstructed signal values at T_TARGET)
    t_sample = (0:length(X_SAMPLE)-1)/F_S;
    xrcon = zeros(size(T_TARGET));
    for k = 1:length(X_SAMPLE)
        xrcon = xrcon + X_SAMPLE(k) * sinc(F_S*(T_TARGET - t_sample(k))); % Add the contribution of the k-th sample to the reconstructed signal (with sinc interpolation)
    end
end
