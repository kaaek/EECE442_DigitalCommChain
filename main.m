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
% % TO-DO: Comment on the figures in the report



% ----------------------------- Fourier -----------------------------
% (a) For a chunked version x(t) choose an appropriate value of T and compute the Fourier series approximation ˆx(t) for different values of n.
% fsConstT([100, 105, 110]); % The optimal option goes last
% (c) Repeat the experiment by varying T while keeping n sufficiently large.
% fsConstN([1.5, 2, 2.2]);
% TO-DO: comment on all the plots.



% ----------------------------- Quantizer -----------------------------
% [t, xt, f_c] = exampleSpeechWave(1);
% % twoLvlQuan(t, xt)
% % % TO-DO: comment on distortion.
% M_values = [4, 32, 128];
% MSE_values = zeros(1, length(M_values));
% [xq, MSE] = uniformQuan(M_values, t, xt);
% TO-DO: report thresholds and Compare the mean square error (MSE) between different quantization rates.

% ----------------------------- Lloyd-Max -----------------------------
M = 20;                             % number of quantization levels
[t, x, f_c] = exampleSpeechWave(1);
[thr, lvl, xq] = lloydMax(x, M);

figure;
plot(t, x, 'r--', 'DisplayName', 'Original Signal');
hold on;
stem(t, xq, 'b-', 'LineStyle', 'None', 'DisplayName', 'Lloyd-Max Quantized');
hold on;
xlabel('time (s)');
ylabel('Amplitude (V)');
legend show;
grid on;
disp(['Mean Squared Error: ', num2str(MSE)]);
