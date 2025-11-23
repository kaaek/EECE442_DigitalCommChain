function ber_mlse = simulate_qpsk_isi_viterbi(Nbits, EbN0dB)

% QPSK ISI Channel + Viterbi Equalizer (MLSE)
% Uses your functions:
%   random_bitstream.m
%   qpsk_mod.m
%   qpsk_demod.m

    if nargin < 1, Nbits = 1e5; end
    if nargin < 2, EbN0dB = 0:2:12; end

    % Ensure even #bits for QPSK
    if mod(Nbits,2) ~= 0
        Nbits = Nbits - 1;
    end

    % ===== Generate random bits using YOUR function =====
    bits_logical = random_bitstream(Nbits);
    bits = double(bits_logical(:));        % column vector

    Nsym = Nbits/2;

    % ===== QPSK modulation using YOUR qpsk_mod =====
    tx_symbols = zeros(Nsym,1);
    for k = 1:Nsym
        b1 = bits(2*k-1);
        b2 = bits(2*k);
        bitpair = sprintf('%d%d', b1, b2);
        tx_symbols(k) = qpsk_mod(bitpair, 1);   % A = 1 because Es=1
    end

    % QPSK constellation from your modulator function:
    % Just generate all 4 symbols:
    const = [
        qpsk_mod('00',1);
        qpsk_mod('01',1);
        qpsk_mod('11',1);
        qpsk_mod('10',1)
    ];

    % ===== Viterbi parameters =====
    M = 4;     % number of states = |S|
    ber_mlse = zeros(1, length(EbN0dB));

    fprintf("\n=== QPSK ISI + Viterbi Using Your Functions ===\n");
    fprintf("Bits = %d | Symbols = %d\n", Nbits, Nsym);
    fprintf("Eb/N0(dB)     BER\n");
    fprintf("---------------------\n");

    Eb = 0.5;   % because Es = 1 and 2 bits/symbol

    % ===== Sweep Eb/N0 =====
    for idx = 1:length(EbN0dB)

        EbN0 = 10^(EbN0dB(idx)/10);
        N0 = Eb / EbN0;
        sigma = sqrt(N0/2);

        % ===== Pass through ISI channel y = 2x_i + x_{i-1} + n_i =====
        y = zeros(Nsym,1);
        prev = 0;   % assume x_0 = 0

        for i = 1:Nsym
            n = sigma*(randn + 1j*randn);
            y(i) = 2*tx_symbols(i) + prev + n;
            prev = tx_symbols(i);
        end

        % ===== Viterbi Equalizer =====

        % Path metrics
        PM = zeros(M,1);
        prev_state = zeros(Nsym, M);

        % ---- Time i=1 ----
        for j = 1:M
            s_prev = 0;         % assumed x_0 = 0
            s_curr = const(j);
            PM(j) = abs(y(1) - (2*s_curr + s_prev))^2;
        end

        % ---- Time i=2..N ----
        for i = 2:Nsym
            PMnew = inf(M,1);
            ps = zeros(M,1);

            for j = 1:M
                s_curr = const(j);
                best_metric = inf;
                best_prev = 1;

                for k = 1:M
                    s_prev = const(k);
                    metric = PM(k) + abs(y(i) - (2*s_curr + s_prev))^2;
                    if metric < best_metric
                        best_metric = metric;
                        best_prev = k;
                    end
                end

                PMnew(j) = best_metric;
                ps(j) = best_prev;
            end

            PM = PMnew;
            prev_state(i,:) = ps;
        end

        % ===== Traceback =====
        [~, st] = min(PM);
        st_seq = zeros(Nsym,1);

        for i = Nsym:-1:1
            st_seq(i) = st;
            st = prev_state(i, st);
            if st == 0, break; end
        end

        x_hat = const(st_seq);

        % ===== Demap using YOUR qpsk_demod =====
        bits_hat = zeros(Nbits,1);

        for k = 1:Nsym
            bstr = qpsk_demod(x_hat(k));  % returns e.g. '10'
            bits_hat(2*k-1) = str2double(bstr(1));
            bits_hat(2*k)   = str2double(bstr(2));
        end

        % ===== BER =====
        ber_mlse(idx) = mean(bits ~= bits_hat);
        fprintf('%8.2f     %.6g\n', EbN0dB(idx), ber_mlse(idx));
    end

    fprintf("====================================\n");

    % ===== Plot =====
    figure;
    semilogy(EbN0dB, ber_mlse, 'o-','LineWidth',1.5);
    grid on;
    xlabel('E_b / N_0 (dB)');
    ylabel('BER');
    title('QPSK ISI Channel with Viterbi Equalizer');
    legend('MLSE');
end
