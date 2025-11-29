function ber_bpsk = simulate_bpsk_awgn(BITSTREAM_LENGTH, SNRdB)
% SIMULATE_BPSK_AWGN 
%   Runs a complete BPSK AWGN simulation using custom mod/demod.
%   Generates random bits, modulates, adds AWGN, demodulates, 
%   computes BER and plots results.
%
%   INPUTS:
%       BITSTREAM_LENGTH : number of bits (default = 1e6)
%       SNRdB            : vector of SNR points in dB (default = 0:2:12)
%
%   OUTPUT:
%       ber_bpsk         : BER values at each SNR

    % ---------- Handle defaults ----------
    if nargin < 1
        BITSTREAM_LENGTH = 1e6;
    end
    if nargin < 2
        SNRdB = 0:2:12;
    end

    % ---------- Generate random bits ----------
    bits_logical = random_bitstream(BITSTREAM_LENGTH);  % 1xBITSTREAM_LENGTH logical array
    bits = double(bits_logical(:));                     % numeric column vector: size BITSTREAM_LENGTHx1
    Nbits = length(bits);

    % ---------- Initialize BER result ----------
    ber_bpsk = zeros(1, length(SNRdB));

    fprintf('\n===== BPSK BER Simulation =====\n');
    fprintf('Bits: %d\n', Nbits);
    fprintf('SNR(dB)        BER\n');
    fprintf('-----------------------------\n');

    % ---------- SNR Loop ----------
    for k = 1:length(SNRdB)
        snr_db = SNRdB(k);

        % Compute Es and amplitude A
        Es = 10^(snr_db/10);
        A  = sqrt(Es);

        % ----- BPSK Modulation (uses your bpsk_mod) -----
        tx = zeros(Nbits,1);
        for n = 1:Nbits
            tx(n) = bpsk_mod(bits(n), A);
        end

        % ----- Add complex AWGN: alpha,beta ~ N(0,0.5) -----
        noise = sqrt(0.5) * (randn(Nbits,1) + 1i*randn(Nbits,1));
        rx = tx + noise;

        % ----- Demodulation (uses your bpsk_demod) -----
        bits_hat = zeros(Nbits,1);
        for n = 1:Nbits
            bits_hat(n) = bpsk_demod(real(rx(n)));
        end

        % ----- BER Calculation -----
        ber_bpsk(k) = mean(bits ~= bits_hat);

        % ----- Labelled Console Output -----
        fprintf('%6d        %.6f\n', snr_db, ber_bpsk(k));
    end

    fprintf('===============================\n\n');

    % ---------- Plot BER ----------
    figure;
    semilogy(SNRdB, ber_bpsk, 'o-', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('BER', 'FontSize', 12);
    title('BPSK BER over AWGN', 'FontSize', 14);
    legend('BPSK', 'Location', 'southwest');

end
