% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function fourier_analysis(SIGNALDURATION, CARRIERFREQUENCY)

fprintf('===============================================\n');
fprintf('   1.2 Fourier Series Experiments running...\n');
fprintf('===============================================\n');

[t, xt, f_max] = exampleSpeechWave(SIGNALDURATION, CARRIERFREQUENCY);
n_opt = ceil(f_max*SIGNALDURATION);
n_1   = ceil(0.9*n_opt);
n_2   = ceil(0.93*n_opt);

% Compute Fourier series approximations for different n values
n_values = [n_1, n_2, n_opt];
fsConstT(t, xt, f_max, n_values, SIGNALDURATION);

t_opt = SIGNALDURATION;
t_1 = 0.4*SIGNALDURATION;
t_2 = 1.2*SIGNALDURATION;

T_values = [t_1, t_opt, t_2];
fsConstN(t, xt, f_max, n_opt, T_values)

% * Discuss trade-off: resolution (n) vs. periodicity assumption (T).
% * Mention effect of complex coefficients; justify taking `real(x̂)` for practical reconstruction.
% * Discuss error energy
end