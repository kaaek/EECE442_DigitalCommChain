% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function coding_decoding_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels)
% CODING_DECODING_ANALYSIS Analyzes full communication chain with Huffman coding
%   Performs quantization (Uniform and Lloyd-Max), Huffman encoding/decoding,
%   and evaluates end-to-end performance metrics including MSE, throughput,
%   and compression efficiency.
%
%   Parameters:
%       SIGNAL_DURATION - Duration of the signal in seconds
%       MESSAGE_FREQUENCY - Frequency of the message signal in Hz
%       CARRIER_FREQUENCY - Frequency of the carrier signal in Hz
%       TARGET_MSE - Target mean squared error for Lloyd-Max quantizer
%       quantizer_levels - Array of quantization levels to test (e.g., [4, 8, 16, 32])

    % Generate signal and sample at Nyquist rate
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
        [xq_u, lvl_u, thr_u, ~] = uniformQuan(t_sample, x_sample, q);
        [thr_lm, lvl_lm, xq_lm, ~] = lloydMax(x_sample, q, TARGET_MSE);
        
        % Compute discrete MSE for quantization (sample-based, not integral-based)
        MSE_quant_u = mean((x_sample(:) - xq_u(:)).^2);
        MSE_quant_lm = mean((x_sample(:) - xq_lm(:)).^2);

        % -------- LOSSLESS CODING (Huffman) --------
        encoded_u = baseline_huffman_V2(xq_u);
        encoded_lm = baseline_huffman_V2(xq_lm);

        % Decode symbols - decoded_symbols is already a string array
        % Convert each string element to a number
        decoded_u = arrayfun(@str2double, encoded_u.decoded_symbols);
        decoded_lm = arrayfun(@str2double, encoded_lm.decoded_symbols);

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
end
