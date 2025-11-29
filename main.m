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
addpath('Analog\PSK\');

[SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels, block_sizes] = constants();

%% Phase 1
% % ==================== 1 - SAMPLING & FOURIER ====================
% sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% fourier_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% % ==================== 2 - QUANTIZATION ====================
% quantization_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% % ==================== 3 - CODING, DECODING... ====================
coding_decoding_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels)

%% Phase 2
% % ==================== 4 - MODULATION DEMO ====================
% modulation_demo()
% simulate_bpsk_awgn()
% simulate_bpsk_coded_awgn()
% simulate_qpsk_awgn()
% simulate_qpsk_coded_awgn()
% simulate_qpsk_isi_viterbi()
% % ==================== 5 - BER ESTIMATION ====================
% BER_estimation_bpsk_awgn()
% BER_estimation_qpsk_isi_mlse()