% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

clc, clearvars;
addpath('util')


% ----------------------------- Sampling  -----------------------------
% f0 = 10;
% example1(f0);
% 
% % TO-DO: Comment on the figures in the report



% ----------------------------- Fourier -----------------------------
% (a) For a chunked version x(t) choose an appropriate value of T and compute the Fourier series approximation ˆx(t) for different values of n.
% fsConstT([100, 105, 110]); % The optimal option goes last
% (c) Repeat the experiment by varying T while keeping n sufficiently large.
% fsConstN([1.5, 2, 2.2]);
% 
% TO-DO: comment on all the plots.



% ----------------------------- Quantizer -----------------------------
% [t, xt, f_c] = exampleSpeechWave(1);
% % twoLvlQuan(t, xt)
% % % TO-DO: comment on distortion.
% M_values = [4, 32, 128];
% MSE_values = zeros(1, length(M_values));
% [xq, MSE] = uniformQuan(M_values, t, xt);
% 
% % TO-DO: report thresholds and Compare the mean square error (MSE) between different quantization rates.
% 
% % ----------------------------- Lloyd-Max -----------------------------
% M = 16;                             % number of quantization levels
% tgtMSE = 0.001;
% [t, x, f_c] = exampleSpeechWave(10);
% [thr, lvl, xq] = lloydMax(x, M, tgtMSE);
% 
% figure;
% plot(t, x, 'r--', 'DisplayName', 'Original Signal');
% hold on;
% stem(t, xq, 'b-', 'LineStyle', 'None', 'DisplayName', 'Lloyd-Max Quantized');
% hold on;
% xlabel('time (s)');
% ylabel('Amplitude (V)');
% legend show;
% grid on;

% TO-DO: * CHANGE STOP CONDITION TO REACHING A MINIMUM MSE
% (g) Integrate your quantization methods into the end-to-end simulation chain.
% (h) Compare the reconstructed signals with different quantization methods.
% (i) Discuss trade-offs between rate (number of levels), distortion (MSE), and source
% statistics.


% ----------------------------- Full Chain ----------------------------- %
% (1) Prepare the input (analog generation):
[t, xt, f_max] = exampleSpeechWave(2, 50); % twenty-one second long speech signal on a "carrier" frequency of 160 Hz (somewhere between a high-pitched male and a low-pitched female voice)
f_Nyquist = 2*f_max; % minimum required for perfect reconstruction.
fs = f_Nyquist;

% (2) Convert from continuous to discrete x[n]:
[t_sample, x_sample] = sample(t, xt, fs);

% (3) Apply quantization to the sampled signal:
M = 16; % number of quantization levels
[xq, MSE] = uniformQuan(M, t_sample, x_sample);

% (4) Apply Lloyd-Max quantization to the sampled signal:
tgtMSE = 0.0000001;
[t_lm, lvl, xq_lm] = lloydMax(x_sample, M, tgtMSE);

% (5) Perform lossless compression using baseline Huffman coding (deterministic step):
C_uniform = baseline_huffman_V2(xq);
C_lloyd   = baseline_huffman_V2(xq_lm);

% (6) Decompress
D_uniform = decode_stream(C_uniform.encoded_bitstring, C_uniform.dict);
D_lloyd = decode_stream(C_lloyd.encoded_bitstring, C_lloyd.dict);
x_hat_uniform = double(D_uniform);
x_hat_lloyd  = double(D_lloyd);

MSE_uniform = mean((xt-x_hat_uniform).^2);
MSE_lloyd = mean((xt-x_hat_lloyd).^2);

% (7) Plot the reconstructed signals for comparison
figure;

% Subplot for Uniform Quantized Signal
subplot(2, 1, 1);
plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
hold on;
plot(t_sample, x_hat_uniform, 'b-', 'DisplayName', 'Uniform Quantized Signal');
xlabel('Time (s)');
ylabel('Amplitude (V)');
legend show;
grid on;
title('Uniform Quantization Reconstruction');

% Subplot for Lloyd-Max Quantized Signal
subplot(2, 1, 2);
plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
hold on;
plot(t_sample, x_hat_lloyd, 'g-', 'DisplayName', 'Lloyd-Max Quantized Signal');
xlabel('Time (s)');
ylabel('Amplitude (V)');
legend show;
grid on;
title('Lloyd-Max Quantization Reconstructed');

% (8) Plot the Mean Square Errors
figure;
bar([MSE_uniform, MSE_lloyd], 'FaceColor', 'flat');
set(gca, 'XTick', 1:2, 'XTickLabel', {'Uniform Quantization', 'Lloyd-Max Quantization'});
ylabel('Mean Square Error (MSE)');
title('MSE Comparison between Quantization Methods');
grid on;
legend('MSE Values', 'Location', 'Best');
