% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function fsConstN(t, xt, f_max, N, T_values)
    % % ---------------------- Init ----------------------
    % duration = 2;                                      % Any number
    % [t, xt, f_c, f_max] = exampleSpeechWave(duration);        % x(t): speech wave
    % E = zeros(1,length(T_values));                     % Init energy array
    % % ------------------ Control Case ------------------
    % T = duration;                                      % T = duration of speech signal to reconstruct the entire clip.
    % N = ceil(f_max * T);                               % Enough harmonics to capture the carrier
    % xhat1 = fs(xt, t, N, T);
    % figure;
    % plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
    % hold on;
    % plot(t, xhat1, 'b-', 'DisplayName', 'FS approx.');
    % title(sprintf('[Example] Control Case for T = duration = %.1f | N = %d', T, N));
    % xlabel('Time (s)');
    % ylabel('Amplitude (V)');
    % legend show;
    % grid on;
    % % --------------------------------------------------
    figure('Name','Fourier Approximation (Fixed Coeffiecients n)');
    for i = 1:length(T_values)
        T_i = T_values(i);
        xhat_i = fs(xt, t, N, T_i);

        E_i = errEn(t, xt, xhat_i); % Needed for the energy plot later.
        E(i) = E_i;
        
        subplot(length(T_values), 1, i);
        plot(t, xt, 'r--', 'DisplayName', 'x(t)', 'LineWidth', 1.5, 'Color', [1 0 0 0.5]);
        hold on;
        plot(t, xhat_i, 'b-', 'DisplayName', sprintf('x\\^(t); T = %.1f', T_i));
        title(sprintf('x\\^(t) (T = %d, N = %.1f)', T_i, N));
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend show;
        grid on;
    end
        
    plotErrEnConstN(t, xt, f_max, N);
    
end

function plotErrEnConstN(t, xt, f_max, N)

    T_values = 0:0.05:ceil(N/f_max);
    for i = 1:length(T_values)
        T_i = T_values(i);
        xhat_i = fs(xt, t, N, T_i);
        E_i = errEn(t, xt, xhat_i);                             % Needed for the energy plot later.
        E(i) = E_i;
    end

    if length(E) > 1
        figure('Name','Error Energy VS Period T');
        plot(T_values, E, 'o-');
        title('Error Energy vs Period T');
        xlabel('T (s)');
        ylabel('Error Energy (E)');
        grid on;
    else
        fprintf('Error Energy: E = %.4f\n', E);
    end
    
end