function x = bpsk_mod(b, A)
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