% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function fsConstN(T_values)
    % ---------------------- Init ----------------------
    duration = 2;                                      % Any number
    [t, xt, f_c, f_max] = exampleSpeechWave(duration);        % x(t): speech wave
    E = zeros(1,length(T_values));                     % Init energy array
    % ------------------ Control Case ------------------
    T = duration;                                      % T = duration of speech signal to reconstruct the entire clip.
    N = ceil(f_max * T);                               % Enough harmonics to capture the carrier
    xhat1 = fs(xt, t, N, T);
    figure;
    plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
    hold on;
    plot(t, xhat1, 'b-', 'DisplayName', 'FS approx.');
    title(sprintf('Control Case for T = duration = %.1f | N = %d', T, N));
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show;
    grid on;
    % --------------------------------------------------
    figure;
    for i = 1:length(T_values)
        T_i = T_values(i);
        xhat_i = fs(xt, t, N, T_i);

        E_i = errEn(t, xt, xhat_i); % Needed for the energy plot later.
        E(i) = E_i;
        
        subplot(length(T_values), 1, i);
        plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
        hold on;
        plot(t, xhat_i, 'b-', 'DisplayName', sprintf('FS approx. for T = %.1f', T_i));
        title(sprintf('Case for T = %.1f | N = %d', T_i, N));
        xlabel('Time (s)');
        ylabel('Amplitude (v)');
        legend show;
        grid on;
    end

    figure;
    plot(T_values, E, 'o-', 'DisplayName', 'Error Energy vs T values');
    title('Error Energy vs T values');
    xlabel('T values');
    ylabel('Error Energy (E)');
    legend show;
    grid on;
end