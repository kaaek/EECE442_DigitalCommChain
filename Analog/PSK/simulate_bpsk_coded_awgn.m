function ber_bpsk_coded = simulate_bpsk_coded_awgn(BITSTREAM_LENGTH, SNRdB)
% SIMULATE_BPSK_CODED_AWGN 
%   BPSK over AWGN with repetition-3 channel coding.
%   - Generates random bits using random_bitstream
%   - Encodes with repetition-3 (0->000, 1->111)
%   - Modulates with your bpsk_mod
%   - Adds complex AWGN (alpha,beta~N(0,0.5))
%   - Demodulates with your bpsk_demod
%   - Decodes by majority vote
%   - Computes BER vs SNR and plots the curve
%
%   INPUTS:
%       BITSTREAM_LENGTH : number of info bits (default = 1e6)
%       SNRdB            : vector of SNR values in dB (default = 0:2:12)
%
%   OUTPUT:
%       ber_bpsk_coded   : coded BER at each SNR

    % ---------- Defaults ----------
    if nargin < 1
        BITSTREAM_LENGTH = 1e6;
    end
    if nargin < 2
        SNRdB = 0:2:12;
    end

    % ---------- Generate original (information) bits ----------
    bits_logical = random_bitstream(BITSTREAM_LENGTH);
    bits = double(bits_logical(:));      % column vector 0/1
    Nbits = length(bits);

    % ---------- Channel encoding: repetition-3 ----------
    % Each bit -> [b b b]
    coded_bits = repelem(bits, 3);       % length = 3*Nbits
    Ncoded = length(coded_bits);

    ber_bpsk_coded = zeros(1, length(SNRdB));

    fprintf('\n===== BPSK Coded (Repetition-3) BER Simulation =====\n');
    fprintf('Info bits: %d   Coded bits: %d\n', Nbits, Ncoded);
    fprintf('SNR(dB)        BER (coded)\n');
    fprintf('--------------------------------\n');

    % ---------- SNR Loop ----------
    for k = 1:length(SNRdB)
        snr_db = SNRdB(k);

        % Es = 10^(SNR/10), A = sqrt(Es)
        Es = 10^(snr_db/10);
        A  = sqrt(Es);

        % ----- BPSK Modulation on coded bits -----
        tx = zeros(Ncoded,1);
        for n = 1:Ncoded
            tx(n) = bpsk_mod(coded_bits(n), A);
        end

        % ----- Add complex AWGN: alpha,beta ~ N(0,0.5) -----
        noise = sqrt(0.5) * (randn(Ncoded,1) + 1i*randn(Ncoded,1));
        rx = tx + noise;

        % ----- BPSK Demodulation -----
        coded_bits_hat = zeros(Ncoded,1);
        for n = 1:Ncoded
            coded_bits_hat(n) = bpsk_demod(real(rx(n)));
        end

        % ----- Channel Decoding: majority vote on groups of 3 -----
        % Reshape into Nbits x 3 matrix
        coded_bits_hat_mat = reshape(coded_bits_hat, 3, []).';
        % Count number of 1s in each triplet
        sum_ones = sum(coded_bits_hat_mat, 2);
        % Majority: >=2 ones -> 1, else 0
        decoded_bits = double(sum_ones >= 2);

        % ----- BER w.r.t. original info bits -----
        ber_bpsk_coded(k) = mean(decoded_bits ~= bits);

        fprintf('%6d        %.6f\n', snr_db, ber_bpsk_coded(k));
    end

    fprintf('================================\n\n');

    % ---------- Plot coded BER ----------
    figure;
    semilogy(SNRdB, ber_bpsk_coded, '^-', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('BER', 'FontSize', 12);
    title('BPSK with Repetition-3 Coding over AWGN', 'FontSize', 14);
    legend('BPSK coded (rep-3)', 'Location', 'southwest');
end
