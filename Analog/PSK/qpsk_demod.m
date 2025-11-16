function b = qpsk_demod(x)
b = zeros(size(x)); % Initialize output array
for k = 1:length(x)
    if real(x(k)) >= 0 && imag(x(k)) >= 0
        b = '11';
    elseif real(x(k)) < 0 && imag(x(k)) >= 0
        b = '10';
    elseif real(x(k)) < 0 && imag(x(k)) < 0
        b = '00';
    elseif real(x(k)) >= 0 && imag(x(k)) < 0
        b = '01';
    end
end
end