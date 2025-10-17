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

SIGNALDURATION = 2;
CARRIERFREQUENCY = 50;

% sampling_analysis(SIGNALDURATION, CARRIERFREQUENCY);
% fourier_analysis(SIGNALDURATION, CARRIERFREQUENCY);
% quantization_analysis(SIGNALDURATION, CARRIERFREQUENCY);

[t, xt, f_max] = exampleSpeechWave(SIGNALDURATION, CARRIERFREQUENCY);
f_Nyquist = 2*f_max;
% ----------------------------------
[t_sample_nqyuist, x_sample_nyquist] = sample(t, xt, f_Nyquist);
% ----------------------------------
M = 64;
tgtMSE = 0.01;
[xq_u, MSEu] = uniformQuan(M, t_sample_nqyuist, x_sample_nyquist, false);
[thr_n, lvl_n, xq_lm, MSElm] = lloydMax(t_sample_nqyuist, x_sample_nyquist, M, tgtMSE);
% ----------------------------------
coding_analysis(xq_u);
coding_analysis(xq_lm);
% ----------------------------------
block_coding_analysis(xq_u , length(unique(xq_u)), [1 2 3 4], true);
block_coding_analysis(xq_lm , length(unique(xq_lm)), [1 2 3 4], true);