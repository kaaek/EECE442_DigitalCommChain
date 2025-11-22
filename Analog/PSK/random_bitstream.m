function bitstream = random_bitstream(BITSTREAM_LENGTH)
arguments
    BITSTREAM_LENGTH (1,1) {mustBePositive, mustBeInteger} = 10e6; % Default length is 10e6 bits
end
random_numbers = rand(1, BITSTREAM_LENGTH); % Generate random numbers between 0 and 1
% Convert to bits by thresholding at 0.5
% Values < 0.5 become 0, values >= 0.5 become 1
bitstream = (random_numbers >= 0.5); 
% fprintf('Random Bitstream:');
% disp(bitstream);
end