function b = bpsk_demod(x)
assert(isnumeric(x) && isscalar(x), 'Input x must be a single numeric value.');
if x >= 0
    b = 1;
elseif x < 0
    b = 0;
end
end