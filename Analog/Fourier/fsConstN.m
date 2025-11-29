function fsConstN(t, xt, f_max, N, T_values)
    % FSCONSTN Plots Fourier series approximation with fixed coefficients N and varying period T
    %   Analyzes how approximation quality changes with period while keeping coefficients fixed
    %   Input: t (time), xt (signal), f_max (max frequency), N (num coefficients), T_values (periods to test)
    figure('Name','Fourier Approximation (Fixed Coefficients N)');
    for i = 1:length(T_values)
        T_i = T_values(i);
        xhat_i = fs(xt, t, N, T_i);        
        subplot(length(T_values), 1, i);
        plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
        hold on;
        plot(t, xhat_i, '-', 'DisplayName', sprintf('Approximation (T=%.1f)', T_i), 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
        title(sprintf('Fourier Series: N = %d, T = %.1f s', N, T_i), 'FontSize', 11);
        xlabel('Time (s)', 'FontSize', 10);
        ylabel('Amplitude', 'FontSize', 10);
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        set(gca, 'FontSize', 9);
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
        title('Error Energy vs. Period', 'FontSize', 12);
        xlabel('T (s)', 'FontSize', 11);
        ylabel('Error Energy', 'FontSize', 11);
        grid on;
        set(gca, 'FontSize', 10);
    else
        fprintf('Error Energy: E = %.4f\n', E);
    end
    
end