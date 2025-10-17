% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function fsConstT(t, xt, f_max, n_values, T)
    % % ---------------------- Init ----------------------
    % duration = 2;
    % [t, xt, f_c, f_max] = exampleSpeechWave(duration);          % x(t): speech wave
    % T = duration;                                               % T = duration of speech signal to reconstruct the entire clip.
    E = zeros(1,length(n_values));                              % Init energy array
    % % ------------------ Optimal Case ------------------
    % n_optimal = ceil(f_max * T);                                % Enough harmonics to capture the entire message's frequency content.
    % xhat1 = fs(xt, t, n_optimal, T);
    % figure;
    % plot(t, xt, 'r--', 'DisplayName', 'x(t)');
    % hold on;
    % plot(t, xhat1, 'b-', 'DisplayName', 'x^(t)');
    % title(sprintf('[Example] Control Case for n = n_{optimal} = %.1f | T = %.1f', n_optimal, T));
    % xlabel('Time (s)');
    % ylabel('Amplitude (V)');
    % legend show;
    % grid on;
    % % --------------------------------------------------
    figure('Name','Fourier Approximation (Fixed Period T)');
    for i = 1:length(n_values)
        n_i = n_values(i);
        xhat_i = fs(xt, t, n_i, T);

        E_i = errEn(t, xt, xhat_i);                             % Needed for the energy plot later.
        E(i) = E_i;
        
        subplot(length(n_values), 1, i);
        plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]);
        hold on;
        plot(t, xhat_i, 'b-', 'DisplayName', sprintf('x\\^(t) for n = %d', n_i));
        title(sprintf('x\\^(t) (n = %d, T = %.1f)', n_i, T));
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend show;
        grid on;
    end

    plotErrEnConstT(t, xt, f_max, T);

end

function plotErrEnConstT(t, xt, f_max, T)

    n_values = 0:1:ceil(f_max*T);
    for i = 1:length(n_values)
        n_i = n_values(i);
        xhat_i = fs(xt, t, n_i, T);
        E_i = errEn(t, xt, xhat_i);                             % Needed for the energy plot later.
        E(i) = E_i;
    end

    if length(E) > 1
        figure('Name','Error Energy VS Coefficients n');
        plot(n_values, E, 'o-');
        title('Error Energy VS Fourier Coefficients n');
        xlabel('n');
        ylabel('Error Energy (E)');
        grid on;
    else
        fprintf('Error Energy: E = %.4f\n', E);
    end
    
end