function [SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels, block_sizes] = constants()
SIGNAL_DURATION   = 2;  
MESSAGE_FREQUENCY = 5;
CARRIER_FREQUENCY = 50;    
TARGET_MSE        = 0.1;  
quantizer_levels  = [4, 8, 16, 32];
block_sizes       = [1 2 3 4];
end