function [lvl, thr] = lloydMaxInit(x_samples, M)
    lvl = linspace(min(x_samples), max(x_samples), M);
    thr = (lvl(1:end-1)+lvl(2:end))/2;          % Define regions according to the selected representation points (with implicit t_0 = -inf and t_M = +inf)
end