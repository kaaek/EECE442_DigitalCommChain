function [T, x_sample1, fs1, x_sample2, fs2] = example1(f0) % Sample a pure (cosine) tone at frequencies 0.5*fN and 2*fN example.
    T = 0:0.001:1;          % Time axis for the original signal (approximately continuous.)
    X = cos(2*pi*f0*T);     % Array of cosine values.
    fN = 2*f0;              % Nyquist's sampling frequency is twice the frequency we care to retain (f0 in this case.)
    fs1 = 0.5*fN;
    fs2 = 2*fN;
    [t_sample1, x_sample1] = sample(T, X, fs1);
    [t_sample2, x_sample2] = sample(T, X, fs2);

    figure;
    subplot(2, 1, 1);
    plot(t_sample1, x_sample1, 'o-', 'DisplayName', sprintf('fs = %.1f Hz', fs1));
    hold on;
    plot(T, X, 'r-', 'DisplayName', 'Original Signal');
    title('Sampling at fs = 0.5 * fN');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show;
    grid on;

    subplot(2, 1, 2);
    plot(t_sample2, x_sample2, 'o-', 'DisplayName', sprintf('fs = %.1f Hz', fs2));
    hold on;
    plot(T, X, 'r--', 'DisplayName', 'Original Signal');
    title('Sampling at fs = 2 * fN');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show;
    grid on;
    % ------------------------ Reconstruction  ------------------------
    [n1, xrcon1] = reconstruct(T, x_sample1, fs1);
    [n2, xrcon2] = reconstruct(T, x_sample2, fs2);

    figure;
    subplot(2, 1, 1);
    plot(T, xrcon1, 'b-', 'DisplayName', 'Reconstructed Signal (fs = 0.5 * fN)');
    hold on;
    stem(n1*(1/fs1), x_sample1, 'r', 'DisplayName', 'Samples');
    title('Reconstruction at fs = 0.5 * fN');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show; grid on;

    subplot(2, 1, 2);
    plot(T, xrcon2, 'b-', 'DisplayName', 'Reconstructed Signal (fs = 2 * fN)');
    hold on;
    stem(n2*(1/fs2), x_sample2, 'r', 'DisplayName', 'Samples');
    title('Reconstruction at fs = 2 * fN');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show; grid on;
end