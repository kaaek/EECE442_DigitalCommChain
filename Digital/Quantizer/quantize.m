% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function xq = quantize(x, thr, lvl)
    xq = zeros(size(x));
    for i = 1:length(x)
        x_i = x(i);
        assigned = false;
        for j = 1:length(thr)
            % Check if the current input value is less than or equal to the current threshold
            if x_i <= thr(j)
                xq(i) = lvl(j);
                assigned = true;
                break;
            end
        end
        % Assign the last level if no thresholds were met
        if ~assigned
            xq(i) = lvl(end);
        end
    end
end