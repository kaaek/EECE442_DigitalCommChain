function ber_qpsk = simulate_qpsk_awgn(BITSTREAM_LENGTH, SNRdB)
% SIMULATE_QPSK_AWGN 
%   Runs a complete QPSK AWGN simulation using your custom qpsk_mod/qpsk_demod.
%   Generates random bits, maps them into QPSK symbols, adds AWGN, 
%   demodulates, computes BER and plots results.
%
%   INPUTS:
%       BITSTREAM_LENGTH : number of bits (default = 1e6)
%       SNRdB            : vector of SNR points in dB (default = 0:2:12)
%
%   OUTPUT:
%       ber_qpsk         : BER values at each SNR

    % ---------- Handle defaults ----------
    if nargin < 1
        BITSTREAM_LENGTH = 1e6;
    end
    if nargin < 2
        SNRdB = 0:2:12;
    end

    % ---------- Generate random bits ----------
    bits_logical = random_bitstream(BITSTREAM_LENGTH);
    bits = double(bits_logical(:));         % numeric column vector 0/1
    Nbits = length(bits);

    % Ensure even number of bits (2 bits per QPSK symbol)
    if mod(Nbits, 2) ~= 0
        bits = bits(1:end-1);
        Nbits = Nbits - 1;
    end
    Nsym = Nbits / 2;

    % ---------- Initialize BER result ----------
    ber_qpsk = zeros(1, length(SNRdB));

    fprintf('\n===== QPSK BER Simulation =====\n');
    fprintf('Bits: %d  (Symbols: %d)\n', Nbits, Nsym);
    fprintf('SNR(dB)        BER\n');
    fprintf('-----------------------------\n');

    % ---------- SNR Loop ----------
    for k = 1:length(SNRdB)
        snr_db = SNRdB(k);

        % Es = 10^(SNR/10), A = sqrt(Es)
        Es = 10^(snr_db/10);
        A  = sqrt(Es);

        % ----- QPSK Modulation using your qpsk_mod -----
        tx = zeros(Nsym,1);
        for n = 1:Nsym
            b1 = bits(2*n-1);           % first bit of symbol
            b2 = bits(2*n);             % second bit of symbol
            bstr = sprintf('%d%d', b1, b2);   % e.g. '01'
            tx(n) = qpsk_mod(bstr, A);       % your QPSK modulator
        end

        % ----- Add complex AWGN: alpha,beta ~ N(0,0.5) -----
        noise = sqrt(0.5) * (randn(Nsym,1) + 1i*randn(Nsym,1));
        rx = tx + noise;

        % ----- QPSK Demodulation using your qpsk_demod -----
        bits_hat = zeros(Nbits,1);
        for n = 1:Nsym
            bstr_hat = qpsk_demod(rx(n));   % e.g. '10'
            % Convert '10' -> [1; 0]
            bits_hat(2*n-1) = str2double(bstr_hat(1));
            bits_hat(2*n)   = str2double(bstr_hat(2));
        end

        % ----- BER Calculation -----
        ber_qpsk(k) = mean(bits ~= bits_hat);

        % ----- Labelled Console Output -----
        fprintf('%6d        %.6f\n', snr_db, ber_qpsk(k));
    end

    fprintf('===============================\n\n');

    % ---------- Plot BER ----------
    figure;
    semilogy(SNRdB, ber_qpsk, 's-', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('BER', 'FontSize', 12);
    title('QPSK BER over AWGN', 'FontSize', 14);
    legend('QPSK', 'Location', 'southwest');

end
