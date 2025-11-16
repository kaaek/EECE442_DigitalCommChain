% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function fourier_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)

fprintf('===============================================\n');
fprintf('   1.2 Fourier Series Experiments running...\n');
fprintf('===============================================\n');

[t, xt, f_max]  = AMWave(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY);
n_opt           = ceil(f_max*SIGNAL_DURATION);
n_1             = ceil(0.9*n_opt);
n_2             = ceil(0.93*n_opt);

% Compute Fourier series approximations for different n values
n_values = [n_1, n_2, n_opt];
fsConstT(t, xt, f_max, n_values, SIGNAL_DURATION);

t_opt = SIGNAL_DURATION;
t_1 = 0.4*SIGNAL_DURATION;
t_2 = 1.2*SIGNAL_DURATION;

T_values = [t_1, t_opt, t_2];
fsConstN(t, xt, f_max, n_opt, T_values)

% * Discuss trade-off: resolution (n) vs. periodicity assumption (T).
% * Mention effect of complex coefficients; justify taking `real(x̂)` for practical reconstruction.
% * Discuss error energy
end