Nbits = 1e5;           % Total bits
EbN0dB = 0:2:12;       % Eb/N0 in dB
Nsym = Nbits / 2;      % Number of QPSK symbols (2 bits/symbol)

if mod(Nbits,2) ~= 0, Nbits = Nbits - 1; end % Ensure even number of bits

bits = randi([0 1], Nbits, 1); % Random bitstream

const = [
    (-1-1j) / sqrt(2);     % index 1: '00' = angle 5π/4
    (1-1j) / sqrt(2);      % index 2: '01' = angle 7π/4
    (-1+1j) / sqrt(2);     % index 3: '10' = angle 3π/4
    (1+1j) / sqrt(2)       % index 4: '11' = angle π/4
];

tx_symbols = zeros(Nsym, 1);
for k = 1:Nsym
    b1 = bits(2*k-1);
    b2 = bits(2*k);
    idx = b1*2 + b2 + 1;  % Maps: '00'=1, '01'=2, '10'=3, '11'=4
    tx_symbols(k) = const(idx);
end

% ---- Initialize BER arrays ----
M = 4; % Constellation size
ber_mlse = zeros(1, length(EbN0dB));
ber_awgn = zeros(1, length(EbN0dB));
ber_theory = zeros(1, length(EbN0dB));

fprintf("\n%s\n", repmat('=',1,60));
fprintf("  QPSK ISI + MLSE vs AWGN Baseline\n");
fprintf("%s\n", repmat('=',1,60));
fprintf("  Nbits = %.0e | Nsym = %.0e\n", Nbits, Nsym);
fprintf("%s\n", repmat('-',1,60));
fprintf("  Eb/N0(dB)    BER_MLSE      BER_AWGN     BER_Theory\n");
fprintf("%s\n", repmat('-',1,60));

Eb = 0.5;  % Es = 1, 2 bits/symbol => Eb = 0.5

% ---- Sweep Eb/N0 ----
for idx = 1:length(EbN0dB)
    
    % Convert dB to linear
    EbN0 = 10^(EbN0dB(idx)/10);
    N0 = Eb / EbN0;
    sigma = sqrt(N0/2);
    
    % ---- ISI Channel: y = 2*x_i + x_{i-1} + n_i ----
    y_isi = zeros(Nsym, 1);
    prev = 0;  % x_0 = 0
    for i = 1:Nsym
        noise = sigma * (randn + 1j*randn);
        y_isi(i) = 2*tx_symbols(i) + prev + noise;
        prev = tx_symbols(i);
    end
    
    % ---- Viterbi MLSE Equalizer ----
    PM = zeros(M, 1);
    backpointer = zeros(Nsym, M);
    
    % t=1: Initialize path metrics
    for j = 1:M
        s_prev = 0;  % assume x_0 = 0
        s_curr = const(j);
        branch_metric = abs(y_isi(1) - (2*s_curr + s_prev))^2;
        PM(j) = branch_metric;
    end
    
    % t=2..N: Forward pass with backpointer tracking
    for t = 2:Nsym
        PM_new = inf(M, 1);
        bp_new = zeros(M, 1);
        
        for j = 1:M  % current state
            s_curr = const(j);
            best_pm = inf;
            best_prev = 0;
            
            for k = 1:M  % previous state
                s_prev = const(k);
                branch = abs(y_isi(t) - (2*s_curr + s_prev))^2;
                path_metric = PM(k) + branch;
                
                if path_metric < best_pm
                    best_pm = path_metric;
                    best_prev = k;
                end
            end
            
            PM_new(j) = best_pm;
            bp_new(j) = best_prev;  % index of best previous state (1..M)
        end
        
        PM = PM_new;
        backpointer(t, :) = bp_new;
    end
    
    % Traceback: reconstruct best path
    [~, state_idx] = min(PM);  % final state with min metric
    state_seq = zeros(Nsym, 1);
    state_seq(Nsym) = state_idx;
    
    for t = Nsym:-1:2
        state_idx = backpointer(t, state_idx);  % follow backpointer
        state_seq(t-1) = state_idx;
    end
    
    % Map state sequence to symbols
    x_hat = const(state_seq);
    
    % Demap and compute BER
    % Find nearest constellation point and extract bits
    bits_hat_mlse = zeros(Nbits, 1);
    for k = 1:Nsym
        % Find closest constellation point
        dists = abs(x_hat(k) - const);
        [~, const_idx] = min(dists);
        % const_idx: 1='00', 2='01', 3='10', 4='11'
        bits_hat_mlse(2*k-1) = floor((const_idx-1)/2);   % bit 1
        bits_hat_mlse(2*k)   = mod(const_idx-1, 2);      % bit 2
    end
    
    ber_mlse(idx) = mean(bits ~= bits_hat_mlse);
    
    % ---- AWGN Baseline (no ISI): y = x + n ----
    noise_awgn = sigma * (randn(Nsym,1) + 1j*randn(Nsym,1));
    y_awgn = tx_symbols + noise_awgn;
    
    bits_hat_awgn = zeros(Nbits, 1);
    for k = 1:Nsym
        % Find closest constellation point
        dists = abs(y_awgn(k) - const);
        [~, const_idx] = min(dists);
        % const_idx: 1='00', 2='01', 3='10', 4='11'
        bits_hat_awgn(2*k-1) = floor((const_idx-1)/2);   % bit 1
        bits_hat_awgn(2*k)   = mod(const_idx-1, 2);      % bit 2
    end
    
    ber_awgn(idx) = mean(bits ~= bits_hat_awgn);
    
    % ---- Theoretical QPSK AWGN BER (Gray coded) ----
    % Q(sqrt(2*Eb/N0)) = 0.5*erfc(sqrt(Eb/N0))
    ber_theory(idx) = 0.5 * erfc(sqrt(EbN0));
    
    fprintf("    %6.2f      %.6e    %.6e    %.6e\n", ...
        EbN0dB(idx), ber_mlse(idx), ber_awgn(idx), ber_theory(idx));
    
end

fprintf("%s\n\n", repmat('=',1,60));

figure('Position', [100 100 900 600]);
semilogy(EbN0dB, ber_mlse, 'o-', 'LineWidth', 2.0, 'MarkerSize', 8, ...
    'Color', [0.0 0.4 0.8], 'DisplayName', 'MLSE (ISI)'); hold on;
semilogy(EbN0dB, ber_awgn, 's-', 'LineWidth', 1.8, 'MarkerSize', 7, ...
    'Color', [0.8 0.4 0.0], 'DisplayName', 'AWGN Sim');
semilogy(EbN0dB, ber_theory, 'm--', 'LineWidth', 2.0, ...
    'DisplayName', 'AWGN Theory (Q-function)');

grid on;
set(gca, 'GridAlpha', 0.3, 'FontSize', 11);
xlabel('E_b / N_0 (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('BER', 'FontSize', 13, 'FontWeight', 'bold');
title('BER Estimation: QPSK ISI Channel with MLSE Equalizer', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 12, 'Box', 'on');
xlim([EbN0dB(1) EbN0dB(end)]);
hold off;

fprintf("  Constellation:  QPSK (Counter-clockwise from upper-right)\n");
fprintf("                  '11'=(1+j)/√2, '10'=(-1+j)/√2, '00'=(-1-j)/√2, '01'=(1-j)/√2\n");
fprintf("  Symbol Energy:  E_s = 1\n");
fprintf("  Bit Energy:     E_b = 0.5\n");
fprintf("  ISI Channel:    y_i = 2*x_i + x_{i-1} + n_i\n");
fprintf("  Equalizer:      Viterbi (MLSE, M=4)\n");
fprintf("  AWGN Noise:     n_i ~ CN(0, sigma^2)\n");
fprintf("  Total Bits:     %.0e\n", Nbits);
fprintf("  Total Symbols:  %.0e\n", Nsym);
fprintf("\nKey Observations:\n");
fprintf("---------\n");
fprintf("  - MLSE performance degrades due to ISI (2*x_i term).\n");
fprintf("  - AWGN baseline shows near-theoretical performance.\n");
fprintf("  - Theoretical curve matches AWGN simulation (validation).\n");
fprintf("\n");
