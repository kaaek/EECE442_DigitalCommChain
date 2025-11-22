% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function fsConstN(t, xt, f_max, N, T_values)
    figure('Name','Fourier Approximation (Fixed Coeffiecients n)');
    for i = 1:length(T_values)
        T_i = T_values(i);
        xhat_i = fs(xt, t, N, T_i);        
        subplot(length(T_values), 1, i);
        plot(t, xt, '--', 'DisplayName', 'x(t)', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
        hold on;
        plot(t, xhat_i, '-', 'DisplayName', sprintf('x^(t); T = %.1f', T_i), 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
        title(sprintf('x^(t) (T = %d, N = %.1f)', T_i, N), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Amplitude (V)', 'FontSize', 11, 'FontWeight', 'bold');
        legend('show', 'Location', 'best', 'FontSize', 9);
        grid on;
        set(gca, 'FontSize', 9, 'LineWidth', 1);
    end
        
    plotErrEnConstN(t, xt, f_max, N);
    
end

function plotErrEnConstN(t, xt, f_max, N)
    T_values = 0:0.05:ceil(N/f_max);
    E = zeros(1,length(T_values)); % Initialize the error energy array to all zeroes
    for i = 1:length(T_values)
        T_i = T_values(i);
        xhat_i = fs(xt, t, N, T_i);
        E_i = errorEnergy(t, xt, xhat_i); % Needed for the energy plot later.
        E(i) = E_i;
    end

    if length(E) > 1
        figure('Name','Error Energy VS Period T');
        plot(T_values, E, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1]);
        title('Error Energy vs Period T', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('T (s)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Error Energy (E)', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca, 'FontSize', 10, 'LineWidth', 1);
    else
        fprintf('Error Energy: E = %.4f\n', E);
    end
    
end