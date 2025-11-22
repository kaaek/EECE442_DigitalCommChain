function x = qpsk_mod(b, A)
if ~(ischar(b) && length(b) == 2 && all(ismember(b, ['0', '1'])))
    error('The string is not a valid two-bit binary string.');
end
switch b
    case '00'
        x = A/sqrt(2) * exp(1i * 5 * pi / 4); % QPSK symbol for '00'
    case '01'
        x = A/sqrt(2) * exp(1i * 7 * pi / 4); % QPSK symbol for '01'
    case '10'
        x = A/sqrt(2) * exp(1i * 3 * pi / 4); % QPSK symbol for '10'
    case '11'
        x = A/sqrt(2) * exp(1i * pi / 4);; % QPSK symbol for '11'
    otherwise
        error('Invalid binary string for QPSK modulation.');
end
end