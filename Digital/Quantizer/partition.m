function lvl = partition(x, thr, M, lvl_prev)
    lvl = zeros(1, M);
    for k = 1:M                                  % For the k-th region
        % Determine the lower bound for the k-th region
        if k == 1
            lower_bound = -inf;  % First region starts from negative infinity
        else
            lower_bound = thr(k-1);  % Lower bound is the previous threshold
        end
        % Determine the upper bound for the k-th region
        if k == M
            upper_bound = inf;  % Last region extends to positive infinity
        else
            upper_bound = thr(k);  % Upper bound is the current threshold
        end
        
        % Extract the samples that fall within the bounds of the k-th region
        region_k = x(x > lower_bound & x <= upper_bound);
        if ~isempty(region_k)                   
            lvl(k) = mean(region_k);            % The level for the k-th region is the mean of the samples in that region
        else
            lvl(k) = lvl_prev(k);
        end
    end
end
