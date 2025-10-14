function xq = twoLvlQuan(t, xt)
    t_1 = mean(xt);
    l_1 = (min(xt)+t_1)/2;
    l_2 = (max(xt)+t_1)/2;
    
    thr = [t_1];
    lvl = [l_1, l_2];

    xq = quan(xt, thr, lvl);
    plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
    hold on;
    stem(t, xq, 'b-', 'LineStyle', 'none', 'DisplayName', 'Quantized Signal');
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    legend show;
    grid on;
end