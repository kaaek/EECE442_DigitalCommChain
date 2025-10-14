function fsConstT(n_values)
    % ---------------------- Init ----------------------
    duration = 2;
    [t, xt, f_c, f_max] = exampleSpeechWave(duration);          % x(t): speech wave
    T = duration;                                               % T = duration of speech signal to reconstruct the entire clip.
    E = zeros(1,length(n_values));                              % Init energy array
    % ------------------ Optimal Case ------------------
    n_optimal = ceil(f_max * T);                                % Enough harmonics to capture the entire message's frequency content.
    xhat1 = fs(xt, t, n_optimal, T);
    figure;
    plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
    hold on;
    plot(t, xhat1, 'b-', 'DisplayName', 'FS approx.');
    title(sprintf('Control Case for n = n_{optimal} = %.1f | T = %.1f', n_optimal, T));
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show;
    grid on;
    % --------------------------------------------------
    figure;
    for i = 1:length(n_values)
        n_i = n_values(i);
        xhat_i = fs(xt, t, n_i, T);
        E_i = errEn(t, xt, xhat_i);                             % Needed for the energy plot later.
        E(i) = E_i;
        subplot(length(n_values), 1, i);
        plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
        hold on;
        plot(t, xhat_i, 'b-', 'DisplayName', sprintf('FS approx. for n = %d', n_i));
        title(sprintf('Case for n = %d | T = %.1f', n_i, T));
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend show;
        grid on;
    end
    figure;
    plot(n_values, E, 'o-', 'DisplayName', 'Error Energy vs n values');
    title('Error Energy vs n values');
    xlabel('n values');
    ylabel('Error Energy (E)');
    legend show;
    grid on;
end