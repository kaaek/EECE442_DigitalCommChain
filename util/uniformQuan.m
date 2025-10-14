function [xq, MSE] = uniformQuan(M, t, xt)
    MSE = zeros(1, length(M));
    figure;
    for i = 1:length(M)
        lvl = linspace(min(xt), max(xt), M(i));
        thr = (lvl(1:end-1) + lvl(2:end)) / 2;
        xq = quan(xt, thr, lvl);
        subplot(length(M), 1, i)
        plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
        hold on;
        stem(t, xq, 'b-', 'LineStyle', 'none', 'DisplayName', 'Quantized Signal');
        title('Original Signal vs Quantized Signal');
        xlabel('Time (s)');
        ylabel('Amplitude (V)');
        legend show;
        grid on;
        MSE(i) = mean((xt - xq).^2);
    end

    figure;
    plot(M,MSE, 'o-', 'DisplayName', 'MSE vs Quantization Levels');
    title('MSE vs Quantization Level M');
    xlabel('M');
    ylabel('MSE');
    legend show;
    grid on;

    % lvl = linspace(min(xt), max(xt), M);
    % thr = (lvl(1:end-1) + lvl(2:end)) / 2;
    % xq = quan(xt, thr, lvl);
    % figure;
    % plot(t, xt, 'r--', 'DisplayName', 'Original Signal');
    % hold on;
    % stem(t, xq, 'b-', 'LineStyle', 'none', 'DisplayName', 'Quantized Signal');
    % title('Original Signal vs Quantized Signal');
    % xlabel('Time (s)');
    % ylabel('Amplitude (V)');
    % legend show;
    % grid on;
    % MSE = mean((xt - xq).^2);
end