function y = sinc(x)
    % SINC Computes sinc function: sin(π*x)/(π*x), with limit sin(0)/0 = 1
    %   Used in signal reconstruction and interpolation via Nyquist sampling theorem
    y = ones(size(x));
    idx = x ~= 0;
    y(idx) = sin(pi*x(idx))./(pi*x(idx));
end