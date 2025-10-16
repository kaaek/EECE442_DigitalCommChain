% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [xq, MSE] = uniformQuan(M, t, xt)
    % uniformQuan quantizes a signal using uniform quantization.
    %
    % Syntax:
    %   [xq, MSE] = uniformQuan(M, t, xt)
    %
    % Inputs:
    %   M   - A vector containing the number of quantization levels.
    %   t   - A vector representing the time instances corresponding to the signal.
    %   xt  - The original signal to be quantized.
    %
    % Outputs:
    %   xq  - The quantized signal.
    %   MSE - A vector containing the Mean Squared Error for each quantization level.
    %
    % Description:
    %   This function performs uniform quantization on the input signal 'xt'
    %   for each value in the vector 'M'. It computes the quantized signal
    %   and the Mean Squared Error (MSE) between the original and quantized
    %   signals. The function also generates plots to visualize the original
    %   and quantized signals, as well as the MSE against the quantization levels.
    MSE = zeros(1, length(M));
    figure;
    for i = 1:length(M) % Plot the quantized signal
        lvl = linspace(min(xt), max(xt), M(i));
        thr = (lvl(1:end-1) + lvl(2:end)) / 2;
        xq = quan(xt, thr, lvl);
        subplot(length(M), 1, i)
        plot(t, xt, 'r--', 'DisplayName', 'x[n]');
        hold on;
        stem(t, xq, 'b-', 'LineStyle', 'none', 'MarkerFaceColor', 'b', 'DisplayName', 'x\^[n]');
        title('[TX] Uniform Quantizer: x[n] vs x\^[n]');
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend show;
        grid on;
        MSE(i) = mean((xt - xq).^2);
    end

    if length(MSE) > 1
        figure;
        plot(M, MSE, 'o-', 'DisplayName', 'MSE vs Quantization Levels');
        title('[TX] Uniform Quantizer: MSE VS Number of Quantization Levels M');
        xlabel('M');
        ylabel('MSE');
        legend show;
        grid on;
    else
        fprintf('[TX] Uniform Quantizer: Mean Squared Error = %.4f\n', MSE);
    end

    % lvl = linspace(min(xt), max(xt), M);
    % thr = (lvl(1:end-1) + lvl(2:end)) / 2;
    % xq = quan(xt, thr, lvl);
    % figure;
    % plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
    % hold on;
    % stem(t, xq, 'b-', 'LineStyle', 'none', 'DisplayName', 'Quantized Signal');
    % title('[TX] Uniform Quantizer: Original Signal VS Quantized Signal');
    % xlabel('Time (s)');
    % ylabel('Amplitude (V)');
    % legend show;
    % grid on;
    % MSE = mean((xt - xq).^2);
end