function [t, xt, f_max] = AMWave(DURATION, F_M, F_C)
    % AMWAVE Generates AM-modulated signal with specified message and carrier frequencies
    %   Computes: xt(t) = (1 + 0.5*m(t)) * cos(2*pi*f_c*t), where m(t) = cos(2*pi*f_m*t)
    %   Output: t (time vector), xt (modulated signal), f_max (maximum frequency component)
    t = 0:0.001:DURATION;
    f_max = F_C + F_M;
    xt = (1 + 0.5*cos(2*pi*F_M*t)) .* cos(2*pi*F_C*t); % Just an AM-modulated wave
end