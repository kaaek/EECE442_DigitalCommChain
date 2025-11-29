function fourier_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)

fprintf(repmat('=',1,70));
fprintf('\n   1.2 Fourier Series Experiments running...\n');
fprintf(repmat('=',1,70));

[t, xt, f_max]  = AMWave(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY);
n_opt           = ceil(f_max*SIGNAL_DURATION);
n_1             = ceil(0.9*n_opt);
n_2             = ceil(0.93*n_opt);
fprintf('\nFourier Series Analysis Results For Constant Signal Duration T = %.2f \n', SIGNAL_DURATION);
fprintf('Fourier Coefficients:\n');
fprintf('n_opt  = %.2f \n', n_opt);
fprintf('n_1    = %.2f \n', n_1);
fprintf('n_2    = %.2f \n', n_2);
n_values = [n_1, n_2, n_opt];
fsConstT(t, xt, f_max, n_values, SIGNAL_DURATION);

t_opt = SIGNAL_DURATION;
t_1 = 0.4*SIGNAL_DURATION;
t_2 = 1.2*SIGNAL_DURATION;
fprintf('Fourier Series Analysis Results For Constant Number of Coefficients N = %.2f \n', n_opt);
fprintf('Durations:\n');
fprintf('t_opt  = %.2f \n', t_opt);
fprintf('t_1    = %.2f \n', t_1);
fprintf('t_2    = %.2f \n', t_2);
T_values = [t_1, t_opt, t_2];
fsConstN(t, xt, f_max, n_opt, T_values)
end