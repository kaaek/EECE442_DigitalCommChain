% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
%
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% * Dark-mode edits by ChatGPT
% ----------------------------------------------------------------------

function coding_analysis(X, OPTS)

fprintf('========================================\n');
fprintf('   3.2 Baseline Huffman is Running...\n');
fprintf('========================================\n');

% ---------- normalize input ----------
if isnumeric(X)
    Xn = string(X(:));
elseif ischar(X) || isstring(X)
    Xn = string(X(:));
else
    error('X must be numeric, char, or string.');
end

if nargin < 2 || ~isstruct(OPTS)
    OPTS = struct();
end

% ---------- run the baseline reference implementation ----------
R = baseline_huffman_V2(Xn, OPTS);

% ---------- DECODE & VERIFY (added) ----------
decoded_msg = R.decoded_symbols;
is_lossless = isequal(string(Xn(:)), decoded_msg);

fprintf('\nDecoded message:\n');
disp(join(decoded_msg,""))   % print as one string (e.g., abfb)
fprintf('Lossless verification (isequal): %d\n', is_lossless);
fprintf('Total coded bits: %d | Bits per symbol: %.4f | Fixed bits/sym: %d\n', ...
        R.total_bits_huffman, R.bits_per_symbol, R.fixed_bits_per_symbol);

% ---------- pull fields for quick access ----------
A_obs   = R.A_observed;
tbl     = R.dict_table;          % Symbol | Probability | Code | Len
H       = R.entropy_bits_per_symbol;
Lfix    = R.fixed_bits_per_symbol;
Lhuff   = R.huffman_avg_bits_per_symbol;
N       = R.N;

% ---------- Plot diagnostics ----------
figure('Name','Baseline Huffman Diagnostics','Color','k'); % black figure

% (1) Probability bar (sorted)
subplot(3,1,1);
[p_sorted, ord] = sort(double(tbl.Probability), 'descend');
sym_sorted      = string(tbl.Symbol(ord));
bar(p_sorted, 'DisplayName','p(a)', 'FaceColor',[0.2 0.6 1]); 
grid on; hold on;
xlabel('Symbol (sorted)', 'Color','w'); 
ylabel('p(a)', 'Color','w');
title('Observed PMF (sorted by probability)', 'FontWeight','bold','Color','w');
xticks(1:numel(sym_sorted));
xticklabels(sym_sorted);
xtickangle(45);
legend show;

% Dark mode styling
set(gca,'Color','k','XColor','w','YColor','w');
legend('TextColor','w','Color','k');

% (2) Code length bar aligned with the same order
subplot(3,1,2);
len_sorted = double(tbl.Len(ord));
bar(len_sorted, 'DisplayName','$\ell(a)$', 'FaceColor',[0.6 0.3 0.9]);  
grid on; hold on;
xlabel('Symbol (sorted)', 'Color','w');
ylabel('bits', 'Color','w');
title('Huffman Code Lengths $\ell(a)$ (aligned with PMF order)', ...
      'Interpreter','latex','FontWeight','bold','Color','w');
xticks(1:numel(sym_sorted));
xticklabels(sym_sorted);
xtickangle(45);
legend('Interpreter','latex','FontSize',9,'Location','best','TextColor','w','Color','k');

% Dark mode styling
set(gca,'Color','k','XColor','w','YColor','w');

% (3) Compare average Huffman length vs. entropy vs. fixed-length
subplot(3,1,3);
hold on; grid on;

% Extract key metrics
H      = R.entropy_bits_per_symbol;
Lbar   = R.huffman_avg_bits_per_symbol;
Lfix   = R.fixed_bits_per_symbol;

% Define values and labels 
vals = [H, Lbar, Lfix];
cats = categorical({'Entropy H(A)', 'Avg Huffman l̄', 'Fixed-length'});
cats = reordercats(cats, {'Entropy H(A)', 'Avg Huffman l̄', 'Fixed-length'});

% Draw bars
bar(cats, vals, 'FaceColor', [0.3 0.6 0.9]);
ylabel('bits / symbol', 'Color','w');
ylim([0, max(vals)*1.35]);
grid on;

% Add numeric labels slightly above bars
for i = 1:numel(vals)
    text(i, vals(i) + 0.04*max(vals), sprintf('%.3f', vals(i)), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontWeight','bold', 'FontSize',10, 'Color','w');
end

title('Comparison of Average Code Lengths: $H(A)$, $\bar{\ell}$, and fixed-length', ...
      'Interpreter','latex', 'FontWeight','bold', 'FontSize',11, 'Color','w');

set(gca, 'FontWeight','bold', 'FontSize',10, 'Color','k', 'XColor','w', 'YColor','w');
legend('TextColor','w','Color','k');
set(gcf,'Color','k'); % keep full figure dark

% ---------- High-level summary --------------------
fprintf('\nBaseline Huffman Report (symbol-wise):\n');
fprintf('--------------------------------------\n');
fprintf('N (number of symbols)       : %d\n', N);
fprintf('|A_observed| / |A_design|   : %d / %d\n', numel(A_obs), numel(R.A_design));
fprintf('Entropy H(A)                : %.4f bits/sym\n', H);
fprintf('Fixed length (design)       : %d bits/sym (total %d)\n', Lfix, R.total_bits_fixed);
fprintf('Huffman avg length          : %.4f bits/sym (total %d)\n', Lhuff, R.total_bits_huffman);
fprintf('Compression vs fixed        : %.2f%%\n', 100*R.compression_gain_vs_fixed);
fprintf('Lossless verified           : %d\n', R.lossless_verified);

% ---------- Explanation for report --------------------
fprintf('\nInterpretation:\n');
fprintf('H(A) represents the theoretical lower bound on bits/symbol.\n');
fprintf('The average Huffman length (l̄) approaches H(A) when symbol probabilities are skewed.\n');
fprintf('When probabilities are uniform, l̄ ≈ fixed-length — meaning no compression benefit.\n');
fprintf('Quantizers with dead-zones or highly uneven histograms yield shorter Huffman codes,\n');
fprintf('because frequent symbols (e.g., zeros) get very short codewords.\n');

end
