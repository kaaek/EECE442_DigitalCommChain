% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
%
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function block_coding_analysis(seq, A, K_values, verify_lossless)

fprintf('========================================\n');
fprintf('   3.2 Block Source Coding is Running...\n');
fprintf('========================================\n');

if nargin < 3 || isempty(K_values), K_values = [1 2 3 4]; end
if nargin < 4, verify_lossless = true; end

% --- Run your original analysis (DO NOT MODIFY original function) ---
results = block_source_coding(seq, A, K_values, verify_lossless);

% --- Show the returned table ---
disp(' ');
disp('==================== BLOCK SOURCE CODING RESULTS ====================');
disp(results);

% --- Establish a K=1 baseline for dependency gains (robust even if K=1 not requested) ---
% Theoretical baseline H1 from symbol-wise entropy; Empirical baseline L1_emp from symbol-wise Huffman
seq_syms = cellstr(string(seq(:)));
R1      = baseline_huffman_V2(seq_syms);       % your existing function
H1      = R1.entropy_bits_per_symbol;          % H(A)
L1_emp  = R1.huffman_avg_bits_per_symbol;      % \bar{\ell} at K=1 (empirical)

% --- Build vectors sorted by K for plotting ---
[Ks, ord]         = sort(results.K);
Hk_per_sym        = results.Hk_per_sym(ord);
L_emp_per_sym     = results.Huff_bits_per_sym_emp(ord);

% --- Dependency gains vs K ---
dep_gain_emp = L1_emp - L_emp_per_sym;   % empirical: \bar{\ell}_{1} - \bar{\ell}_{K}
dep_gain_th  = H1     - Hk_per_sym;      % theoretical: H_1 - H_K/K

% Numerical hygiene (clip tiny negatives from float error)
tol = 1e-12;
dep_gain_emp(abs(dep_gain_emp) < tol) = 0;
dep_gain_th(abs(dep_gain_th)   < tol) = 0;

% ===================== PLOT: Short-range Dependency Gains vs K =====================
figure('Name','Short-range Dependency Gains','Color','w');
plot(Ks, dep_gain_emp, '-o', 'LineWidth', 1.6, 'DisplayName', '$\bar{\ell}_{1}-\bar{\ell}_{K}$'); hold on;
plot(Ks, dep_gain_th,  '--s','LineWidth', 1.6, 'DisplayName', '$H_{1}-H_{K}/K$');
yline(0,'k:','HandleVisibility','off');
grid on;
xlabel('Block size K');
ylabel('bits / symbol');
title('Gains from Short-Range Dependencies vs K', 'FontWeight','bold');
legend('Interpreter','latex','Location','best');

% Mark the best empirical K (>1), if any
mask = Ks > 1;
if any(mask)
    [bestGain, idxLocal] = max(dep_gain_emp(mask));
    idxs  = find(mask);
    kBest = Ks(idxs(idxLocal));
    yBest = bestGain;
    plot(kBest, yBest, 'd', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName','Max empirical gain');
    text(kBest, yBest + 0.03*max([dep_gain_emp; dep_gain_th; 1e-6]), ...
        sprintf('K=%d, +%.3f b/sym', kBest, yBest), ...
        'HorizontalAlignment','center','FontWeight','bold');
end

% ===================== Printed Discussion =====================
fprintf('\n==================== DISCUSSION ====================\n');
for i = 1:numel(Ks)
    fprintf('K=%d | H_K/K=%.4f | l̄_emp=%.4f | DepGain(emp)=+%.4f | DepGain(theory)=+%.4f\n', ...
        Ks(i), Hk_per_sym(i), L_emp_per_sym(i), dep_gain_emp(i), dep_gain_th(i));
end

fprintf(['\n• Positive gains for K>1 indicate **short-range dependencies** (quantizer memory or\n' ...
         '  correlated signal statistics) that make certain blocks more probable, allowing\n' ...
         '  shorter Huffman codewords.\n' ...
         '• When $H_K/K \\approx \\bar{\\ell}_K$, coding is near-optimal for the empirical model.\n' ...
         '• If gains plateau/shrink as K increases, likely causes are sparse block histograms,\n' ...
         '  limited sample size (noisy p estimates), or little correlation beyond short ranges.\n' ...
         '• Dead-zones / repetitive regions often yield larger gains since block PMFs are skewed.\n']);

end
