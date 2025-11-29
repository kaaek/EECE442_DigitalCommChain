function coding_analysis(X, OPTS)

fprintf('========================================\n');
fprintf('   3.2 Baseline Huffman is Running...\n');
fprintf('========================================\n');

% ---------- normalize input ----------
if isnumeric(X)
    Xn = string(X(:));
elseif ischar(X)
    Xn = string(X(:));              % one symbol per character
elseif isstring(X)
    if isscalar(X)
        % split scalar string into characters
        Xn = string(cellstr(char(X)'));  % one symbol per character
    else
        Xn = string(X(:));
    end
else
    error('X must be numeric, char, or string.');
end

% ---------- ensure OPTS exists ----------
if nargin < 2 || ~isstruct(OPTS)
    OPTS = struct();
end

% ---------- run the baseline reference implementation ----------
R = baseline_huffman_V2(Xn, OPTS);

% ---------- DECODE & VERIFY (robust) ----------
decoded_msg = [];
cand_dec = ["decoded_symbols","decoded_sequence","decoded_syms","decoded"];
for f = cand_dec
    if isfield(R, f)
        tmp = R.(f);
        if isstring(tmp) || ischar(tmp)
            decoded_msg = string(tmp(:)); break
        elseif iscell(tmp)
            decoded_msg = string(tmp(:)); break
        elseif isnumeric(tmp)
            decoded_msg = string(tmp(:)); break
        end
    end
end

if isempty(decoded_msg)
    warning('Baseline did not expose decoded output; skipping equality check.');
    is_lossless = NaN;
else
    is_lossless = isequal(string(Xn(:)), decoded_msg(:));
end

% ---------- Bits-per-symbol ------------------
bps = [];
cand_bps = ["bits_per_symbol","huffman_avg_bits_per_symbol"];
for f = cand_bps
    if isfield(R, f), bps = R.(f); break; end
end
if isempty(bps) && isfield(R, "total_bits_huffman")
    bps = double(R.total_bits_huffman) / numel(Xn);
end
if isempty(bps), bps = NaN; end

fprintf('\nDecoded message:\n');
if isempty(decoded_msg)
    disp('(not available from baseline)')
else
    disp(join(decoded_msg,""))   % print as one string
end
fprintf('Lossless verification (isequal): %s\n', string(is_lossless));
fprintf('Total coded bits: %g | Bits per symbol: %.4f | Fixed bits/sym: %g\n', ...
        gf(R,'total_bits_huffman',NaN), bps, gf(R,'fixed_bits_per_symbol',NaN));

% ---------- pull fields for quick access ----------
A_obs   = gf(R,'A_observed',[]);
tbl     = gf(R,'dict_table',table());
H       = gf(R,'entropy_bits_per_symbol',NaN);
Lfix    = gf(R,'fixed_bits_per_symbol',NaN);
Lhuff   = gf(R,'huffman_avg_bits_per_symbol',bps);
N       = gf(R,'N',numel(Xn));

% ---------- Plot diagnostics ----------
figure('Name','Baseline Huffman Diagnostics','Color','k'); % black figure

% (1) Probability bar (sorted)
subplot(3,1,1);
if ~isempty(tbl) && any(strcmpi(tbl.Properties.VariableNames,'Probability'))
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
    lg = legend('show'); if ~isempty(lg), set(lg,'TextColor','w','Color','k'); end
else
    text(0.5,0.5,'No dict\_table.Probability available','Color','w','HorizontalAlignment','center');
end
set(gca,'Color','k','XColor','w','YColor','w');

% (2) Code length bar aligned with the same order
subplot(3,1,2);
if ~isempty(tbl) && any(strcmpi(tbl.Properties.VariableNames,'Len'))
    if ~exist('ord','var'), ord = 1:height(tbl); end
    len_sorted = double(tbl.Len(ord));
    bar(len_sorted, 'DisplayName','$\ell(a)$', 'FaceColor',[0.6 0.3 0.9]);  
    grid on; hold on;
    xlabel('Symbol (sorted)', 'Color','w');
    ylabel('bits', 'Color','w');
    title('Huffman Code Lengths $\ell(a)$ (aligned with PMF order)', ...
          'Interpreter','latex','FontWeight','bold','Color','w');
    if exist('sym_sorted','var')
        xticks(1:numel(sym_sorted));
        xticklabels(sym_sorted);
        xtickangle(45);
    end
    lg = legend('Interpreter','latex','FontSize',9,'Location','best'); 
    if ~isempty(lg), set(lg,'TextColor','w','Color','k'); end
else
    text(0.5,0.5,'No dict\_table.Len available','Color','w','HorizontalAlignment','center');
end
set(gca,'Color','k','XColor','w','YColor','w');

% (3) Compare average Huffman length vs. entropy vs. fixed-length
subplot(3,1,3);
hold on; grid on;

vals = [H, Lhuff, Lfix];
cats = categorical({'Entropy H(A)', 'Avg Huffman l̄', 'Fixed-length'});
cats = reordercats(cats, {'Entropy H(A)', 'Avg Huffman l̄', 'Fixed-length'});

bar(cats, vals, 'FaceColor', [0.3 0.6 0.9]);
ylabel('bits / symbol', 'Color','w');
ylim_hi = max(vals(~isnan(vals)));
if isempty(ylim_hi), ylim_hi = 1; end
ylim([0, ylim_hi*1.35 + eps]);

for i = 1:numel(vals)
    if ~isnan(vals(i))
        text(i, vals(i) + 0.04*ylim_hi, sprintf('%.3f', vals(i)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontWeight','bold', 'FontSize',10, 'Color','w');
    end
end
title('Comparison of Average Code Lengths: $H(A)$, $\bar{\ell}$, and fixed-length', ...
      'Interpreter','latex', 'FontWeight','bold', 'FontSize',11, 'Color','w');

set(gca, 'FontWeight','bold', 'FontSize',10, 'Color','k', 'XColor','w', 'YColor','w');
legend('off');
set(gcf,'Color','k');

% ---------- High-level summary --------------------
fprintf('\nBaseline Huffman Report (symbol-wise):\n');
fprintf('--------------------------------------\n');
fprintf('N (number of symbols)       : %d\n', N);
if isfield(R,'A_design')
    fprintf('|A_observed| / |A_design|   : %d / %d\n', numel(A_obs), numel(R.A_design));
else
    fprintf('|A_observed|                 : %d (A_design not provided)\n', numel(A_obs));
end
fprintf('Entropy H(A)                : %.4f bits/sym\n', H);
fprintf('Fixed length (design)       : %g bits/sym (total %g)\n', Lfix, gf(R,'total_bits_fixed',NaN));
fprintf('Huffman avg length          : %.4f bits/sym (total %g)\n', Lhuff, gf(R,'total_bits_huffman',NaN));
if isfield(R,'compression_gain_vs_fixed')
    fprintf('Compression vs fixed        : %.2f%%\n', 100*R.compression_gain_vs_fixed);
end
if ~isnan(is_lossless)
    fprintf('Lossless verified           : %d\n', is_lossless);
else
    fprintf('Lossless verified           : (not available)\n');
end

% ---------- Explanation for report --------------------
fprintf('\nInterpretation:\n');
fprintf('H(A) represents the theoretical lower bound on bits/symbol.\n');
fprintf('The average Huffman length (l̄) approaches H(A) when symbol probabilities are skewed.\n');
fprintf('When probabilities are uniform, l̄ ≈ fixed-length — meaning no compression benefit.\n');
fprintf('Quantizers with dead-zones or highly uneven histograms yield shorter Huffman codes,\n');
fprintf('because frequent symbols (e.g., zeros) get very short codewords.\n');

end

% ============== local helpers ==============
function v = gf(S, field, defaultVal)
if isstruct(S) && isfield(S, field)
    v = S.(field);
else
    v = defaultVal;
end
end
