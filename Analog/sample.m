% ----------------------------------------------------------------------
% author: Khalil El Kaaki
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors if exists.
% ----------------------------------------------------------------------

% function [tSample, xSample] = sample(xt, t, F)
%     T = 1/F;                                        % Define step size
%     tSample = t(1):T:t(end);                        % Traverse from the signal's support in step size T    
%     xSample = interp1(t, xt, tSample, 'spline');    % 'spline' option found by trial and error. Check interp1 documentation.
% end

function [tSample, xSample] = sample(xt, t, F)
    % xt, t are column or row vectors of same length
    T = 1/F;
    tSample = t(1):T:t(end);

    % Choose a relative tolerance (adjust as needed for your data scale)
    relTol = 1e-6;                       % example relative tolerance
    tol = relTol * max(1, range(t));     % absolute tolerance based on scale

    % Find unique t within tolerance and map indices
    [tUniq, ~, ic] = uniquetol(t(:), tol);

    % Aggregate xt per group (mean). Ensure xt is a column.
    xtCol = xt(:);
    ytAgg = accumarray(ic, xtCol, [], @mean);

    % Sort unique t and aggregated y so tUniq is strictly increasing
    [tSorted, order] = sort(tUniq);
    ySorted = ytAgg(order);

    % Interpolate using the cleaned coordinates
    xSample = interp1(tSorted, ySorted, tSample, 'spline');
end
