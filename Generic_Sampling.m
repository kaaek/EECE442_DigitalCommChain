% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

clc, clearvars



% ----------------------------- Sampling  -----------------------------
% f0 = 10;
% example1(f0);
% TO-DO: Comment on the figures in the report



% ----------------------------- Fourier -----------------------------
% (a) For a chunked version x(t) choose an appropriate value of T and compute the Fourier series approximation ˆx(t) for different values of n.
% fsConstT([90, 95, 100, 125, 300]);
% (c) Repeat the experiment by varying T while keeping n sufficiently large.
% fsConstN([1.5, 1.75, 2, 2.25, 3]);
% TO-DO: comment on all the plots.



% ----------------------------- Quantizer -----------------------------
% [t, xt, f_c] = exampleSpeechWave(1);
% twoLvlQuan(t, xt)
% TO-DO: comment on distortion.



% ----------------------------- Functions -----------------------------
function [t_sample, x_sample] = sample(t, xt, fs)
    %SAMPLE   Sample a signal at a given frequency
    %   [T_SAMPLE, X_SAMPLE] = SAMPLE(T, XT, FS) resamples the input signal
    %   XT defined on the time axis T at a sampling frequency FS.
    %
    %   Inputs:
    %       T   - row vector, time axis
    %       XT  - row vector, signal values
    %       FS  - scalar, sampling frequency
    %
    %   Outputs:
    %       T_SAMPLE - row vector, sampled time axis
    %           By definition, this is: t = n*T_s
    %       X_SAMPLE - row vector, sampled signal values.
    %           By definition this is: x[n] = x(t) at t = nT_s such that n is in the integer set Z.
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       [ts, xs] = sample(t, xt, 20);
    %       plot(ts, xs, 'o')
    Ts = 1/fs;
    t_sample = t(1):Ts:t(end);
    x_sample = interp1(t, xt, t_sample, 'spline'); % 'spline' option found by trial and error. Check interp1 documentation.
end

function y = sinc(x)
    y = ones(size(x));
    idx = x ~= 0;
    y(idx) = sin(pi*x(idx))./(pi*x(idx));
end

function [n, xrcon] = reconstruct(t, x_sample, fs)
    %RECONSTRUCT   Reconstruct a signal from its samples
    %   XRCON = RECONSTRUCT(T, X_SAMPLE, FS) reconstructs the original signal
    %   from the sampled signal X_SAMPLE defined on the time axis T at a 
    %   sampling frequency FS.
    %
    %   Inputs:
    %       T        - row vector, time axis for the original signal
    %       X_SAMPLE - row vector, sampled signal values
    %       FS       - scalar, sampling frequency used for sampling
    %
    %   Output:
    %       XRCON    - row vector, reconstructed signal values
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       fs = 20;
    %       [ts, xs] = sample(t, xt, fs);
    %       xr = reconstruct(t, xs, fs);
    %       plot(t, xt, 'r-', t, xr, 'b--')
    %       legend('Original Signal', 'Reconstructed Signal');
    Ts = 1/fs;
    n = 0:length(x_sample)-1;
    t_sample = n*Ts;              
    xrcon = zeros(size(t));
    for k = 1:length(n)
        xrcon = xrcon + x_sample(k) * sinc(fs * (t - t_sample(k)));
    end
end

function xhat = fs(xt, t, n, T)
    % FS   Fourier Series approximation
    %   [XHAT, CK] = FS(XT, T, N, T) computes the Fourier Series
    %   approximation of the signal XT defined on the time axis T using
    %   N coefficients.
    %
    %   Inputs:
    %       XT - row vector, signal values
    %       T  - row vector, time axis
    %       N  - scalar, number of coefficients to use
    %       T  - scalar, period of the signal
    %
    %   Outputs:
    %       XHAT - row vector, approximated signal values
    %       CK   - row vector, Fourier coefficients
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       n = 5;
    %       T = 1/10; % Period of the signal
    %       [xhat, ck] = fs(xt, t, n, T);
    %       plot(t, xt, 'r-', t, xhat, 'b--')
    %       legend('Original Signal', 'Fourier Series Approximation');
    k = -n:n;
    ck = zeros(1, length(k));
    xhat = zeros(size(t));

    for idx = 1:length(k)   % Calculate coefficients c_k
        ck(idx) = fcoef(t, xt, T, k(idx));
    end
    for idx = 1:length(k)   % Reconstruct signal using c_k
        xhat = xhat + ck(idx) * exp(1i*2*pi*k(idx)*t/T);
    end
end

function cj = fcoef(t, xt, T, j) % Helper function that calculates the jth fourier series coefficient.
    %FCOEF   Fourier series coefficient calculation
    %   CJ = FCOEF(T, XT, T, J) computes the j-th Fourier series coefficient
    %   for the signal XT defined on the time axis T with period T.
    %
    %   Inputs:
    %       T  - row vector, time axis
    %       XT - row vector, signal values
    %       T  - scalar, period of the signal
    %       J  - scalar, index of the coefficient to compute
    %
    %   Output:
    %       CJ - scalar, j-th Fourier series coefficient
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       T = 1/10; % Period of the signal
    %       j = 1;
    %       cj = fcoef(t, xt, T, j);
    %       disp(cj);
    cj = (1/T) * trapz(t, xt .* exp(-1i*2*pi*j*t/T));
end

function E = errEn (t, xt, xhat)
    % ERREREN   Error Energy Calculation
    %   E = ERREREN(T, XT, XHAT) computes the error energy between the 
    %   original signal XT and the approximated signal XHAT over the time 
    %   axis T.
    %
    %   Inputs:
    %       T    - row vector, time axis
    %       XT   - row vector, original signal values
    %       XHAT - row vector, approximated signal values
    %
    %   Output:
    %       E - scalar, computed error energy
    %
    %   Example:
    %       t = 0:0.001:1;
    %       xt = cos(2*pi*10*t);
    %       xhat = xt + 0.1*randn(size(xt)); % Example of an approximated signal
    %       E = errEn(t, xt, xhat);
    integrand = (xt - xhat).^2;
    E = trapz(t, integrand);
end

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
    plot(T, X, 'r-', 'DisplayName', 'Original Signal');
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

function [t, xt, f_c] = exampleSpeechWave(duration)
    t = 0:0.0001:duration;
    f_c = 50;  % Carrier frequency
    f_m = 5;    % Message frequency
    xt = (1 + 0.5*cos(2*pi*f_m*t)) .* cos(2*pi*f_c*t);
end

function fsConstT(n_values)
    % ---------------------- Init ----------------------
    duration = 2;
    [t, xt, f_c] = exampleSpeechWave(duration);        % x(t): speech wave
    T = duration;                                      % T = duration of speech signal to reconstruct the entire clip.
    E = zeros(1,length(n_values));                     % Init energy array
    % ------------------ Control Case ------------------
    f_max = f_c;                                % Max frequency in signal
    n_optimal = ceil(f_max * T);                % Enough harmonics to capture the carrier
    xhat1 = fs(xt, t, n_optimal, T);
    % ---------------------- Figure --------------------
    figure;
    plot(t, xt, 'r-', 'DisplayName', 'Original Signal');
    hold on;
    plot(t, xhat1, 'b--', 'DisplayName', 'FS approx.');
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
        E_i = errEn(t, xt, xhat_i); % Needed for the energy plot later.
        E(i) = E_i;
        subplot(length(n_values), 1, i);
        plot(t, xt, 'r-', 'DisplayName', 'Original Signal');
        hold on;
        plot(t, xhat_i, 'b--', 'DisplayName', sprintf('FS approx. for n = %d', n_i));
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

function fsConstN(T_values)
    % ---------------------- Init ----------------------
    duration = 2;                                      % Any number
    [t, xt, f_c] = exampleSpeechWave(duration);        % x(t): speech wave
    E = zeros(1,length(T_values));                       % Init energy array
    % ------------------ Control Case ------------------
    T = duration;                                      % T = duration of speech signal to reconstruct the entire clip.
    f_max = f_c;                                       % Max frequency in signal
    N = ceil(f_max * T);                               % Enough harmonics to capture the carrier
    xhat1 = fs(xt, t, N, T);
    % ---------------------- Figure --------------------
    figure;
    plot(t, xt, 'r-', 'DisplayName', 'Original Signal');
    hold on;
    plot(t, xhat1, 'b--', 'DisplayName', 'FS approx.');
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
        plot(t, xt, 'r-', 'DisplayName', 'Original Signal');
        hold on;
        plot(t, xhat_i, 'b--', 'DisplayName', sprintf('FS approx. for T = %.1f', T_i));
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

function xq = quan(x, thr, lvl)
    % quan - Quantizes the input signal based on specified thresholds and levels.
    %
    % Syntax:
    %   xq = quan(x, thr, lvl)
    %
    % Inputs:
    %   x   - Input signal to be quantized (vector).
    %   thr - Thresholds for quantization (vector). The length of thr should be one less than the length of lvl.
    %   lvl - Levels corresponding to the thresholds (vector). The length of lvl should be one more than the length of thr.
    %
    % Outputs:
    %   xq  - Quantized output signal (vector) of the same size as input x.
    %
    % Description:
    %   The function quantizes the input signal x based on the provided thresholds and levels. 
    %   Each value in x is replaced by the corresponding level based on the defined thresholds.
    %   If a value in x is less than or equal to the first threshold, it is assigned the first level.
    %   If it is greater than the last threshold, it is assigned the last level. 
    %   Values between thresholds are assigned the corresponding levels.
    xq = zeros(size(x));
    for i = 1:length(x)
        x_i = x(i);
        if x_i <= thr(1)
            xq(i) = lvl(1);
        elseif x_i > thr(end)
            xq(i) = lvl(end);
        else
            for j = 1:length(thr)-1
                if x_i > thr(j) && x_i <= thr(j+1)
                    xq(i) = lvl(j+1);
                    break;
                end
            end
        end
    end
end

function xq = twoLvlQuan(t, xt)
    t_1 = mean(xt);
    l_1 = (min(xt)+t_1)/2;
    l_2 = (max(xt)+t_1)/2;
    
    thr = [t_1];
    lvl = [l_1, l_2];

    xq = quan(xt, thr, lvl);
    plot(t, xt, 'r-', 'DisplayName', 'Original Signal');
    hold on;
    stem(t, xq, 'b--', 'LineStyle', 'none', 'DisplayName', 'Quantized Signal');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show;
    grid on;
end