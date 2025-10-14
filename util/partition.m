% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function lvl = partition(x, thr, M, lvl_prev)
    lvl = zeros(1, M);
    for k = 1:                                  % For the k-th region:
        if k == 1
            lower_bound = -inf;
        else
            lower_bound = thr(k-1);
        end
        if k == M
            upper_bound = inf;
        else
            upper_bound = thr(k);
        end
        region_k = x(x > lower_bound & x <= upper_bound);        
        if ~isempty(region_k)                   
            lvl(k) = mean(region_k);            % The level for the k-th region is the mean of the samples in that region
        else
            lvl(k) = lvl_prev(k);
        end
    end
end
