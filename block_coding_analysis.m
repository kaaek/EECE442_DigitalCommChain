% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
%
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function block_coding_analysis(seq, A, K_values, verify_lossless)
% BLOCK_CODING_ANALYSIS Analyze block source coding performance
%
% This function performs block source coding analysis on a given sequence
% using specified parameters. It computes the empirical and theoretical 
% performance metrics, including entropy and average bits per symbol, 
% and visualizes the results through plots.
%
% Inputs:
%   seq          - Input sequence to be analyzed (vector or matrix).
%   A            - Alphabet size used for coding (scalar).
%   K_values     - Vector of block sizes to analyze (optional).
%                 Default is [1 2 3 4] if not provided.
%   verify_lossless - Boolean flag to verify lossless coding (optional).
%                     Default is true if not provided.
%
% Outputs:
%   The function does not return any outputs but displays results in the 
%   command window and generates plots for visual analysis.
%
% Example:
%   block_coding_analysis(my_sequence, 256, [1 2 3 4], true);
%
% See also: BLOCK_SOURCE_CODING, BASELINE_HUFFMAN_V2

fprintf('========================================\n');
fprintf('   3.2 Block Source Coding is Running...\n');
fprintf('========================================\n');

if nargin < 3 || isempty(K_values), K_values = [1 2 3 4]; end
if nargin < 4, verify_lossless = true; end


results = block_source_coding(seq, A, K_values, verify_lossless);

% --- Show the returned table ---
disp(' ');
disp('==================== BLOCK SOURCE CODING RESULTS ====================');
disp(results);


% Theoretical baseline H1 from symbol-wise entropy; Empirical baseline L1_emp from symbol-wise Huffman
seq_syms = cellstr(string(seq(:)));
R1      = baseline_huffman_V2(seq_syms);       
H1      = R1.entropy_bits_per_symbol;          % H(A)
L1_emp  = R1.huffman_avg_bits_per_symbol;      %  at K=1 (empirical)

% --- Build vectors sorted by K for plotting ---
[Ks, ord]         = sort(results.K);
Hk_per_sym        = results.Hk_per_sym(ord);
L_emp_per_sym     = results.Huff_bits_per_sym_emp(ord);

% --- Dependency gains vs K ---
dep_gain_emp = L1_emp - L_emp_per_sym;  
dep_gain_th  = H1     - Hk_per_sym;      

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
