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

%% Phase 1
% % ==================== 1 - SAMPLING & FOURIER ====================
% sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% fourier_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% % ==================== 2 - QUANTIZATION ====================
% quantization_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY) % Signal generation is implicit in this module
% % ==================== 3 - CODING, DECODING... ====================
% coding_decoding_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels)
%% Phase 2
