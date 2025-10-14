function [t, xt, f_c, f_max] = exampleSpeechWave(duration)
    t = 0:0.0001:duration;
    f_c = 50;  % Carrier frequency
    f_m = 5;    % Message frequency
    f_max = f_c + f_m;
    xt = (1 + 0.5*cos(2*pi*f_m*t)) .* cos(2*pi*f_c*t);
end