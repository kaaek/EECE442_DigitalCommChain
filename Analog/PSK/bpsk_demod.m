function b = bpsk_demod(x)
    % BPSK_DEMOD Demodulates BPSK symbol to binary bit using threshold at zero
    %   Input: x (complex symbol). Output: b (0 if real(x) < 0, else 1)
    assert(isnumeric(x) && isscalar(x), 'Input x must be a single numeric value.');
if x >= 0
    b = 1;
elseif x < 0
    b = 0;
end
end