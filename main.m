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

[SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels, block_sizes] = constants();

% ==================== 1 - SAMPLING & FOURIER ====================
% sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% fourier_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module

% ==================== 2 - QUANTIZATION ====================
% quantization_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module

% ==================== 3 - CODING, DECODING... ====================
[t, xt, f_max] = AMWave(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY);
f_Nyquist = 2*f_max;
[t_sample, x_sample] = sample(t, xt, f_Nyquist);

% Preallocate results tables
num_quantizers = length(quantizer_levels);
results_table_lm = table('Size', [num_quantizers, 6], ...
    'VariableTypes', {'double','double','double','double','double','double'}, ...
    'VariableNames', {'M','MSE_quant','MSE_final','TotalBits','BitsPerSymbol','Latency'});
results_table_u = table('Size', [num_quantizers, 6], ...
    'VariableTypes', {'double','double','double','double','double','double'}, ...
    'VariableNames', {'M','MSE_quant','MSE_final','TotalBits','BitsPerSymbol','Latency'});

idx = 0;
for q = quantizer_levels
    idx = idx + 1;
    % -------- QUANTIZATION --------
    [xq_u, lvl_u, thr_u, MSE_quant_u] = uniformQuan(t_sample, x_sample, q);
    [thr_lm, lvl_lm, xq_lm, MSE_quant_lm] = lloydMax(x_sample, q, TARGET_MSE);

    % -------- LOSSLESS CODING (Huffman) --------
    encoded_u = baseline_huffman_V2(xq_u);
    encoded_lm = baseline_huffman_V2(xq_lm);

    % Decode symbols - decoded_symbols is already a string array
    % Convert each string element to a number
    decoded_u = arrayfun(@str2double, encoded_u.decoded_symbols);
    decoded_lm = arrayfun(@str2double, encoded_lm.decoded_symbols);
    
    % Verify lossless encoding for this iteration
    if idx == 1
        fprintf('Debug first iteration (q=%d):\n', q);
        fprintf('  xq_u size: %dx%d\n', size(xq_u));
        fprintf('  decoded_u size: %dx%d\n', size(decoded_u));
        fprintf('  Original xq_u (first 10): \n'); disp(xq_u(1:min(10,end)));
        fprintf('  Decoded (first 10): \n'); disp(decoded_u(1:min(10,end)));
        fprintf('  Match: %d\n', isequal(xq_u(:), decoded_u(:)));
    end

    % -------- DEQUANTIZATION --------
    % The decoded values should already be the quantization levels
    x_dequant_u = decoded_u(:);  % Ensure column vector
    x_dequant_lm = decoded_lm(:);  % Ensure column vector

    % -------- METRICS --------
    total_bits_u = length(encoded_u.encoded_bitstring);
    bits_per_symbol_u = total_bits_u / length(xq_u);
    latency_symbols_u = 1;
    MSE_final_u = mean((x_sample(:) - x_dequant_u).^2);

    total_bits_lm      = length(encoded_lm.encoded_bitstring);
    bits_per_symbol_lm = total_bits_lm / length(xq_lm);
    latency_symbols_lm = 1;
    MSE_final_lm       = mean((x_sample(:) - x_dequant_lm).^2);

    % Store results in preallocated tables
    results_table_u(idx,:) = {q, MSE_quant_u, MSE_final_u, total_bits_u, bits_per_symbol_u, latency_symbols_u};
    results_table_lm(idx,:) = {q, MSE_quant_lm, MSE_final_lm, total_bits_lm, bits_per_symbol_lm, latency_symbols_lm};
end

% Verify lossless compression
fprintf('\n================= VERIFICATION =================\n');
fprintf('Lossless Compression Check:\n');
fprintf('  All symbols should match after encode/decode cycle\n');
fprintf('  If MSE_quant ≈ MSE_final, the chain is working correctly\n\n');

disp('\n================= FULL CHAIN PERFORMANCE (LLOYD-MAX QUANTIZER) =================\n');
disp(results_table_lm);

disp('\n================= FULL CHAIN PERFORMANCE (UNIFORM QUANTIZER) =================\n');
disp(results_table_u);

% Plot results AFTER the loop completes
figure('Name','MSE VS M (Lloyd-Max)');
plot(results_table_lm.M, results_table_lm.MSE_quant, '-o', 'LineWidth', 2, 'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1], 'DisplayName', 'Quantization MSE'); 
hold on;
plot(results_table_lm.M, results_table_lm.MSE_final, '-s', 'LineWidth', 2, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], 'DisplayName', 'End-to-End MSE');
grid on;
xlabel('Quantizer Levels (M)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Mean Squared Error', 'FontSize', 12, 'FontWeight', 'bold');
title('Lloyd-Max Quantizer + Huffman Coding', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name','Throughput VS M (Lloyd-Max)');
plot(results_table_lm.M, results_table_lm.BitsPerSymbol, '-s', 'LineWidth', 2, 'Color', [0.47 0.67 0.19], 'MarkerFaceColor', [0.47 0.67 0.19]); 
grid on;
xlabel('Quantizer Levels (M)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Bits per Symbol', 'FontSize', 12, 'FontWeight', 'bold');
title('Lloyd-Max Compression Rate', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name','MSE VS M (Uniform)');
plot(results_table_u.M, results_table_u.MSE_quant, '-o', 'LineWidth', 2, 'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1], 'DisplayName', 'Quantization MSE'); 
hold on;
plot(results_table_u.M, results_table_u.MSE_final, '-s', 'LineWidth', 2, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], 'DisplayName', 'End-to-End MSE');
grid on;
xlabel('Quantizer Levels (M)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Mean Squared Error', 'FontSize', 12, 'FontWeight', 'bold');
title('Uniform Quantizer + Huffman Coding', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10, 'LineWidth', 1);

figure('Name','Throughput VS M (Uniform)');
plot(results_table_u.M, results_table_u.BitsPerSymbol, '-s', 'LineWidth', 2, 'Color', [0.47 0.67 0.19], 'MarkerFaceColor', [0.47 0.67 0.19]); 
grid on;
xlabel('Quantizer Levels (M)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Bits per Symbol', 'FontSize', 12, 'FontWeight', 'bold');
title('Uniform Quantizer Compression Rate', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 10, 'LineWidth', 1);