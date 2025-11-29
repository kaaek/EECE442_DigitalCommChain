function b = qpsk_demod(x)
    % QPSK_DEMOD Demodulates complex QPSK symbol to 2-bit binary string using quadrant detection
    %   Maps: Real>=0,Imag>=0->'11'; Real<0,Imag>=0->'10'; Real<0,Imag<0->'00'; Real>=0,Imag<0->'01'
    %   Input: x (complex QPSK symbols). Output: b (2-bit binary strings)
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