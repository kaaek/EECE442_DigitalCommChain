% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------


function coding_analysis(X, OPTS)


fprintf('========================================\n');
fprintf('   3.2 Baseline Huffman is Running...\n');
fprintf('========================================\n');

% ---------- --------normalize input ---------------------------
if isnumeric(X)
    Xn = string(X(:));
elseif ischar(X) || isstring(X)
    Xn = string(X(:));
else
    error('X must be numeric, char, or string.');
end

if nargin < 2
    OPTS = struct();
end

% ---------- run the baseline reference implementation ----------
R = baseline_huffman_V2(Xn, OPTS);

% ---------- pull fields for quick access ----------
A_obs   = R.A_observed;
tbl     = R.dict_table;          % Symbol | Probability | Code | Len
H       = R.entropy_bits_per_symbol;
Lfix    = R.fixed_bits_per_symbol;
Lhuff   = R.huffman_avg_bits_per_symbol;
N       = R.N;

% ---------- Plot diagnostics (3 stacked like your sampling figure) ----------
figure('Name','Baseline Huffman Diagnostics','Color','w');

% (1) Probability bar (sorted)
subplot(3,1,1);
[p_sorted, ord] = sort(tbl.Probability, 'descend');
sym_sorted      = string(tbl.Symbol(ord));
bar(p_sorted, 'DisplayName','p(a)'); grid on; hold on;
xlabel('Symbol (sorted)'); ylabel('p(a)');
title('Observed PMF (sorted by probability)');
xticks(1:numel(sym_sorted));
xticklabels(sym_sorted);
xtickangle(45);
legend show;

% (2) Code length bar aligned with the same order
subplot(3,1,2);
len_sorted = tbl.Len(ord);
bar(len_sorted, 'DisplayName','\ell(a)'); grid on; hold on;
xlabel('Symbol (sorted)'); ylabel('bits');
title('Huffman Code Lengths \ell(a) (aligned with PMF order)');
xticks(1:numel(sym_sorted));
xticklabels(sym_sorted);
xtickangle(45);
legend show;

% (3) Scatter p vs length with Shannon ideal overlay -log2 p
subplot(3,1,3);
scatter(p_sorted, len_sorted, 36, 'filled', 'DisplayName','Huffman \ell(a)'); hold on; grid on;
ideal = -log2(max(p_sorted, eps));
plot(p_sorted, ideal, 'k--', 'DisplayName','Shannon bound -log_2 p(a)', 'LineWidth',1.2);
xlabel('p(a)'); ylabel('bits');
title('Code Optimality Check: \ell(a) vs -log_2 p(a)');
legend show;

% ---------- High-level summary (like your MSE print) ----------
fprintf('\nBaseline Huffman Report (symbol-wise):\n');
fprintf('--------------------------------------\n');
fprintf('N (number of symbols)       : %d\n', N);
fprintf('|A_observed| / |A_design|   : %d / %d\n', numel(A_obs), numel(R.A_design));
fprintf('Entropy H(A)                : %.4f bits/sym\n', H);
fprintf('Fixed length (design)       : %d bits/sym (total %d)\n', Lfix, R.total_bits_fixed);
fprintf('Huffman avg length          : %.4f bits/sym (total %d)\n', Lhuff, R.total_bits_huffman);
fprintf('Compression vs fixed        : %.2f%%\n', 100*R.compression_gain_vs_fixed);
fprintf('Lossless verified           : %d\n', R.lossless_verified);


end


