function ber_qpsk_coded = simulate_qpsk_coded_awgn(BITSTREAM_LENGTH, SNRdB)
% SIMULATE_QPSK_CODED_AWGN 
%   QPSK over AWGN with repetition-3 channel coding.
%   - Generates random bits using random_bitstream
%   - Encodes with repetition-3 (0->000, 1->111)
%   - Maps coded bits to QPSK symbols using qpsk_mod
%   - Adds complex AWGN (alpha,beta~N(0,0.5))
%   - Demodulates using qpsk_demod
%   - Decodes by majority vote
%   - Computes BER vs SNR and plots the curve
%
%   INPUTS:
%       BITSTREAM_LENGTH : number of info bits (default = 1e6)
%       SNRdB            : vector of SNR points in dB (default = 0:2:12)
%
%   OUTPUT:
%       ber_qpsk_coded   : coded BER at each SNR

    % ---------- Defaults ----------
    if nargin < 1
        BITSTREAM_LENGTH = 1e6;
    end
    if nargin < 2
        SNRdB = 0:2:12;
    end

    % ---------- Generate original (information) bits ----------
    bits_logical = random_bitstream(BITSTREAM_LENGTH);
    bits = double(bits_logical(:));         % column vector 0/1
    Nbits = length(bits);

    % Make Nbits a multiple of 6 (LCM of 2 for QPSK, 3 for repetition-3)
    rem6 = mod(Nbits, 6);
    if rem6 ~= 0
        bits = bits(1:end-rem6);
        Nbits = Nbits - rem6;
    end

    % ---------- Channel encoding: repetition-3 ----------
    % Each bit -> [b b b]
    coded_bits = repelem(bits, 3);          % length = 3*Nbits
    Ncoded = length(coded_bits);

    % Sanity: Ncoded is even (2 bits per QPSK symbol) and multiple of 3 (by construction)
    Nsym = Ncoded / 2;                      % QPSK symbols

    ber_qpsk_coded = zeros(1, length(SNRdB));

    fprintf('\n===== QPSK Coded (Repetition-3) BER Simulation =====\n');
    fprintf('Info bits: %d   Coded bits: %d   Symbols: %d\n', Nbits, Ncoded, Nsym);
    fprintf('SNR(dB)        BER (coded)\n');
    fprintf('--------------------------------\n');

    % ---------- SNR Loop ----------
    for k = 1:length(SNRdB)
        snr_db = SNRdB(k);

        % Es = 10^(SNR/10), A = sqrt(Es)
        Es = 10^(snr_db/10);
        A  = sqrt(Es);

        % ----- QPSK Modulation on coded bits -----
        tx = zeros(Nsym,1);
        for n = 1:Nsym
            b1 = coded_bits(2*n-1);          % first bit of symbol
            b2 = coded_bits(2*n);            % second bit of symbol
            bstr = sprintf('%d%d', b1, b2);  % e.g. '01'
            tx(n) = qpsk_mod(bstr, A);       % your QPSK modulator
        end

        % ----- Add complex AWGN: alpha,beta ~ N(0,0.5) -----
        noise = sqrt(0.5) * (randn(Nsym,1) + 1i*randn(Nsym,1));
        rx = tx + noise;

        % ----- QPSK Demodulation -----
        coded_bits_hat = zeros(Ncoded,1);
        for n = 1:Nsym
            bstr_hat = qpsk_demod(rx(n));    % e.g. '10'
            coded_bits_hat(2*n-1) = str2double(bstr_hat(1));
            coded_bits_hat(2*n)   = str2double(bstr_hat(2));
        end

        % ----- Channel Decoding: majority vote on groups of 3 -----
        coded_bits_hat_mat = reshape(coded_bits_hat, 3, []).';   % Nbits x 3
        sum_ones = sum(coded_bits_hat_mat, 2);
        decoded_bits = double(sum_ones >= 2);                    % majority rule

        % ----- BER w.r.t. original info bits -----
        ber_qpsk_coded(k) = mean(decoded_bits ~= bits);

        fprintf('%6d        %.6f\n', snr_db, ber_qpsk_coded(k));
    end

    fprintf('================================\n\n');

    % ---------- Plot coded BER ----------
    figure;
    semilogy(SNRdB, ber_qpsk_coded, 'd-', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('BER', 'FontSize', 12);
    title('QPSK with Repetition-3 Coding over AWGN', 'FontSize', 14);
    legend('QPSK coded (rep-3)', 'Location', 'southwest');

end
