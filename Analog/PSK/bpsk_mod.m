function x = bpsk_mod(b, A)
    % BPSK_MOD Modulates binary bit (0 or 1) to BPSK constellation point
    %   Input: b (0 or 1), A (amplitude). Output: x (complex BPSK symbol: -A for b=0, +A for b=1)
if isempty(b)
    error('Input binary sequence b cannot be empty.');
end

if b == 0
    x = exp(1i * pi);
elseif b == 1
    x = 1;
else
    error('Input binary sequence b must be either 0 or 1.');
end
x = x * A;
end