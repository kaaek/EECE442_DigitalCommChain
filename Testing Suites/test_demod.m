
const = [
    (1+1j)/sqrt(2) * exp(1j*5*pi/4);  % 1: '00'
    (1-1j)/sqrt(2) * exp(1j*7*pi/4);  % 2: '01'
    (-1+1j)/sqrt(2) * exp(1j*3*pi/4); % 3: '10'
    (-1-1j)/sqrt(2) * exp(1j*pi/4)    % 4: '11'
];

fprintf("Constellation test:\n");
for i = 1:4
    fprintf("  const(%d) = %.4f + %.4fj\n", i, real(const(i)), imag(const(i)));
end

% Test bit-to-symbol mapping
bits = [0 0 0 1 1 0 1 1];
fprintf("\nBit-to-symbol mapping:\n");
for k = 1:4
    b1 = bits(2*k-1);
    b2 = bits(2*k);
    idx = b1*2 + b2 + 1;
    fprintf("  bits[%d:%d] = [%d %d] => idx=%d => symbol %.4f+%.4fj\n", ...
        2*k-1, 2*k, b1, b2, idx, real(const(idx)), imag(const(idx)));
end

% Test demodulation (nearest neighbor)
fprintf("\nDemodulation (nearest neighbor):\n");
test_symbols = const + 0.1*(randn(4,1) + 1j*randn(4,1));
for k = 1:4
    dists = abs(test_symbols(k) - const);
    [~, const_idx] = min(dists);
    b1 = floor((const_idx-1)/2);
    b2 = mod(const_idx-1, 2);
    fprintf("  symbol %d: closest idx=%d => [%d %d]\n", k, const_idx, b1, b2);
end

fprintf("\nTest completed successfully.\n");
