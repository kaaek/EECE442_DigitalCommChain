function fsConstT(t, xt, f_max, n_values, T)
    figure('Name','Fourier Approximation (Fixed Period T)');
    for i = 1:length(n_values) % one test per n value.
        n_i = n_values(i);
        xhat_i = fs(xt, t, n_i, T);
        subplot(length(n_values), 1, i);
        plot(t, xt, '--', 'DisplayName', 'Original', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]);
        hold on;
        plot(t, xhat_i, '-', 'DisplayName', sprintf('Approximation (n=%d)', n_i), 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
        title(sprintf('Fourier Series: n = %d, T = %.1f s', n_i, T), 'FontSize', 11);
        xlabel('Time (s)', 'FontSize', 10);
        ylabel('Amplitude', 'FontSize', 10);
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        set(gca, 'FontSize', 9);
    end
    nSweep(t, xt, f_max, T); % Plotting the error energy vs many values of coefficients n.
end

function nSweep(t, xt, f_max, T)
    n_values = 0:1:ceil(f_max*T);
    E = zeros(1,length(n_values)); % Initialize the error energy array to all zeroes
    for i = 1:length(n_values)
        n_i = n_values(i);
        xhat_i = fs(xt, t, n_i, T);
        E_i = errorEnergy(t, xt, xhat_i); % Needed for the energy plot later.
        E(i) = E_i;
    end

    if length(E) > 1
        figure('Name','Error Energy VS Coefficients n');
        plot(n_values, E, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1]);
        title('Error Energy vs. Number of Coefficients', 'FontSize', 12);
        xlabel('n', 'FontSize', 11);
        ylabel('Error Energy', 'FontSize', 11);
        grid on;
        set(gca, 'FontSize', 10);
    else
        fprintf('Error Energy: E = %.4f\n', E);
    end
    
end