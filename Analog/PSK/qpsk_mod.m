function x = qpsk_mod(b, A)
    % QPSK_MOD Modulates 2-bit binary string to QPSK constellation point with unit energy
    %   Counter-clockwise from upper-right: '11'→(1+j), '10'→(-1+j), '00'→(-1-j), '01'→(1-j)
    %   Input: b (2-bit string '00','01','10','11'), A (amplitude). Output: x (QPSK symbol)
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
        x = A/sqrt(2) * exp(1i * pi / 4); % QPSK symbol for '11'
    otherwise
        error('Invalid binary string for QPSK modulation.');
end
end