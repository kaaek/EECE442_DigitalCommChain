% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

clc; clearvars;

% Add to path the directories containing scripts (do not change!)
addpath('Testing Suites\')
addpath('Metrics\');
addpath('Math\');
addpath('Digital\Quantizer');
addpath('Digital\Quantizer\Uniform Quantizer\');
addpath('Digital\Quantizer\Lloyd-Max Quantizer\');
addpath('Digital\Encoder\');
addpath('Analog\');
addpath('Analog\Fourier\');

[SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels, block_sizes] = constants()

% ==================== 1 - SAMPLING & FOURIER ====================

sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module

fourier_analysis(SIGNAL_DURATION, CARRIER_FREQUENCY) % Signal generation is implicit in this module

% ==================== 2 - QUANTIZATION ====================

quantization_analysis(SIGNAL_DURATION, CARRIER_FREQUENCY) % Signal generation is implicit in this module

% % ==================== 3 - CODING, DECODING... ====================
[t, xt, f_max] = exampleSpeechWave(SIGNAL_DURATION, CARRIER_FREQUENCY);
f_Nyquist = 2*f_max;
[t_sample, x_sample] = sample(t, xt, f_Nyquist);
results_table_lm = table();
results_table_u  = table();

for q = quantizer_levels
    % -------- QUANTIZATION --------
    [thr_u, lvl_u, xq_u, MSE] = uniformQuan(q, x_sample, t_sample, false);
    [thr_lm, lvl_lm, xq_lm, MSE_quant] = lloydMax(x_sample, q, TARGET_MSE);

    % -------- LOSSLESS CODING (Huffman) --------
    encoded_u = baseline_huffman_V2(xq_u);
    encoded_lm = baseline_huffman_V2(xq_lm);

    % Decode symbols
    decoded_symbols_u = encoded_u.decoded_symbols;
    [~, decoded_idx_u] = ismember(decoded_symbols_u, string(lvl_u));

    decoded_symbols_lm = encoded_lm.decoded_symbols;
    [~, decoded_idx_lm] = ismember(decoded_symbols_lm, string(lvl_lm));

    % -------- DEQUANTIZATION --------
    x_dequant_u = lvl_u(decoded_idx_u);
    x_dequant_lm = lvl_lm(decoded_idx_lm);

    % -------- METRICS --------
    total_bits_u = length(encoded_u.encoded_bitstring);
    bits_per_symbol_u = total_bits_u / length(xq_u);
    latency_symbols_u = 1;
    MSE_final_u = mean((x_sample - x_dequant_lm).^2);

    total_bits_lm      = length(encoded_lm.encoded_bitstring);
    bits_per_symbol_lm = total_bits_lm / length(xq_lm);
    latency_symbols_lm = 1;
    MSE_final_lm       = mean((x_sample - x_dequant_lm).^2);

    results_table_u = [results_table_u;
        table(q, MSE, MSE_final_u, total_bits_u, bits_per_symbol_u, latency_symbols_u, ...
        'VariableNames',{'M','MSE_quant','MSE_final','TotalBits','BitsPerSymbol','Latency'})];

    results_table_lm = [results_table_lm; 
        table(q, MSE_quant, MSE_final_lm, total_bits_lm, bits_per_symbol_lm, latency_symbols_lm, ...
        'VariableNames', {'M','MSE_quant','MSE_final','TotalBits','BitsPerSymbol','Latency'})];
end

disp('================= FULL CHAIN PERFORMANCE (LLOYD-MAX QUANTIZER) =================');
disp(results_table_lm);

disp('================= FULL CHAIN PERFORMANCE (UNIFORM QUANTIZER) =================');
disp(results_table_u);

figure('Name','MSE VS M (Lloyd-Max)');
plot(results_table_lm.M, results_table_lm.MSE_final, '-o','LineWidth',1.5); grid on;
xlabel('Quantizer Levels (M)'); ylabel('End-to-End MSE');
title('Quantizer + Huffman Performance (Lloyd-Max)');

figure('Name','Throughput VS M (Lloyd-Max)');
plot(results_table_lm.M, results_table_lm.BitsPerSymbol, '-s','LineWidth',1.5); grid on;
xlabel('Quantizer Levels (M)'); ylabel('Bits per Symbol');
title('Throughput (bits/symbol) per Quantizer (Lloyd-Max)');

figure('Name','MSE VS M (Uniform)');
plot(results_table_u.M, results_table_u.MSE_final, '-o','LineWidth',1.5); grid on;
xlabel('Quantizer Levels (M)'); ylabel('End-to-End MSE');
title('Quantizer + Huffman Performance (Uniform)');

figure('Name','Throughput VS M (Uniform)');
plot(results_table_u.M, results_table_u.BitsPerSymbol, '-s','LineWidth',1.5); grid on;
xlabel('Quantizer Levels (M)'); ylabel('Bits per Symbol');
title('Throughput (bits/symbol) per Quantizer (Uniform)');