function [xq, lvl, thr, MSEE] = uniformQuan(t, xt, M)    
    lvl = linspace(min(xt), max(xt), M);
    thr = (lvl(1:end-1) + lvl(2:end)) / 2; % Calculate the quantization levels and thresholds for uniform quantization
    xq = quantize(xt, thr, lvl); % Quantize the original signal using the calculated thresholds and levels
    MSEE = MSE(xt, xq, t);
    
    % fprintf('\nUNIFORM QUANTIZER (M = %s)', num2str(M));
    % fprintf('\nQUANTIZATION LEVELS:\n');
    % fprintf('%.4f ', lvl);
    % fprintf('\nTHRESHOLDS:\n');
    % fprintf('%.4f ', thr);
    % fprintf('\nMEAN SQUARED ERROR (MSE): %.4e\n\n', MSEE);
end