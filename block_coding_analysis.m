% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
%
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function [results_uniform, results_lloydmax, compare_table] = block_coding_analysis( ...
    seq_uniform, seq_lloydmax, A, K_values, verify_lossless)

fprintf('========================================\n');
fprintf('   3.2 Block Source Coding is Running...\n');
fprintf('========================================\n');

if nargin < 4 || isempty(K_values), K_values = [1 2 3 4]; end
if nargin < 5, verify_lossless = true; end
if isempty(seq_lloydmax)
    seq_lloydmax = seq_uniform;
    fprintf('Note: Lloyd–Max sequence not provided; using Uniform sequence for both.\n');
end

% Run analyses
results_uniform  = block_source_coding(seq_uniform,  A, K_values, verify_lossless);
results_lloydmax = block_source_coding(seq_lloydmax, A, K_values, verify_lossless);

% Latency metric
results_uniform.Latency_symbols  = results_uniform.K;
results_lloydmax.Latency_symbols = results_lloydmax.K;

% Cross-comparison
[Ks_common, iu, il] = intersect(results_uniform.K, results_lloydmax.K, 'stable');
compare_table = {};
for t = 1:numel(Ks_common)
    K = Ks_common(t);
    compare_table(end+1,:) = {"Uniform",  K, results_uniform.Total_bits(iu(t)), ...
        results_uniform.Bits_per_symbol(iu(t)), results_uniform.Latency_symbols(iu(t))}; %#ok<AGROW>
    compare_table(end+1,:) = {"LloydMax", K, results_lloydmax.Total_bits(il(t)), ...
        results_lloydmax.Bits_per_symbol(il(t)), results_lloydmax.Latency_symbols(il(t))}; %#ok<AGROW>
end

% Print tables
disp(' ');
disp('==================== UNIFORM RESULTS ====================');
disp(results_uniform);
disp('==================== LLOYD–MAX RESULTS ==================');
disp(results_lloydmax);

disp('==================== UNIFORM vs LLOYD–MAX COMPARISON ====================');
cmp_hdr = {'Quantizer','K','Total_bits','Bits_per_symbol','Latency_symbols'};
cmp_tbl = cell2table(compare_table, 'VariableNames', cmp_hdr);
disp(cmp_tbl);

% Baselines for Uniform (ensure string array for baseline_huffman_V2)
seq_syms_u = string(seq_uniform(:));
R1_u = baseline_huffman_V2(seq_syms_u);
H1_u = R1_u.entropy_bits_per_symbol;
L1_emp_u = R1_u.huffman_avg_bits_per_symbol;

% Vectors for plots (Uniform path)
[Ks_u, ord_u]   = sort(results_uniform.K);
Hk_per_sym_u    = results_uniform.Hk_per_sym(ord_u);
L_emp_per_sym_u = results_uniform.Huff_bits_per_sym_emp(ord_u);
dep_gain_emp_u  = L1_emp_u - L_emp_per_sym_u;
dep_gain_th_u   = H1_u     - Hk_per_sym_u;

% ================== PLOT 1: Bits per symbol vs K (dark mode) ==================
f1 = figure('Name','Bits per symbol vs K (Uniform vs Lloyd–Max)','Color','k');
ax1 = axes('Parent',f1); set(ax1,'Color','k','XColor','w','YColor','w'); hold on; grid on; ax1.GridColor=[.5 .5 .5];

[Ks_plot_u, ord_plot_u] = sort(results_uniform.K);
bps_u = results_uniform.Bits_per_symbol(ord_plot_u);
plot(Ks_plot_u, bps_u, '-o','LineWidth',1.6,'Color',[0 0.6 1],'DisplayName','Uniform');

[Ks_plot_l, ord_plot_l] = sort(results_lloydmax.K);
bps_l = results_lloydmax.Bits_per_symbol(ord_plot_l);
plot(Ks_plot_l, bps_l, '--s','LineWidth',1.6,'Color',[1 0.3 0.3],'DisplayName','Lloyd–Max');

xlabel('Block size K','Color','w'); ylabel('bits / symbol','Color','w');
title('Throughput (bits per symbol) vs K — Uniform vs Lloyd–Max','Color','w','FontWeight','bold');
legend('TextColor','w','Color',[0.1 0.1 0.1]);

% Inline labels
text(Ks_plot_u(end), bps_u(end), '  Uniform',   'Color','w','FontWeight','bold','VerticalAlignment','middle');
text(Ks_plot_l(end), bps_l(end), '  Lloyd–Max','Color','w','FontWeight','bold','VerticalAlignment','middle');

% ============ PLOT 2: Dependency gains vs K (dark mode, Uniform) =============
f2 = figure('Name','Short-range Dependency Gains (Uniform)','Color','k');
ax2 = axes('Parent',f2); set(ax2,'Color','k','XColor','w','YColor','w'); hold on; grid on; ax2.GridColor=[.5 .5 .5];

plot(Ks_u, dep_gain_emp_u, '-o','LineWidth',1.6,'Color',[0 0.6 1],'DisplayName','$\bar{\ell}_{1}-\bar{\ell}_{K}$');
plot(Ks_u, dep_gain_th_u,  '--s','LineWidth',1.6,'Color',[1 0.8 0.2],'DisplayName','$H_{1}-H_{K}/K$');
yline(0,'Color',[0.7 0.7 0.7],'LineStyle',':','HandleVisibility','off');
xlabel('Block size K','Color','w'); ylabel('bits / symbol','Color','w');
title('Gains from Short-Range Dependencies (Uniform)','Color','w','FontWeight','bold');
legend('TextColor','w','Color',[0.1 0.1 0.1],'Interpreter','latex','Location','best');
text(Ks_u(end), dep_gain_emp_u(end), '  Uniform','Color','w','FontWeight','bold','VerticalAlignment','middle');

% ===================== Printed summary =====================
fprintf('\n==================== SUMMARY: Uniform vs Lloyd–Max ====================\n');
nrows = size(compare_table, 1);
for t = 1:2:(nrows-1)
    Ku = compare_table{t,2};
    bpsU = compare_table{t,4};
    bpsL = compare_table{t+1,4};
    fprintf('K=%d | Uniform=%.4f b/sym | Lloyd–Max=%.4f b/sym | Δ=%.4f b/sym\n', Ku, bpsU, bpsL, bpsU-bpsL);
end

% ===================== HUFFMAN TREES (dark mode, labeled) =====================
% Choose a K to visualize (use last K that both results computed)
if ~isempty(K_values)
    K_show = K_values(end);
    try
        local_plot_huffman_tree_for(seq_uniform,  A, K_show, 'Uniform');
        local_plot_huffman_tree_for(seq_lloydmax, A, K_show, 'Lloyd–Max');
    catch ME
        warning('Huffman tree plotting failed for K=%d: %s', K_show, ME.message);
    end
end

end % ===== end main function =====


% ========================= LOCAL HELPERS =========================

function local_plot_huffman_tree_for(seq, A, K, labelName)
    % Build K-block symbols
    s = seq(:)';  N = numel(s);
    numBlocks = floor(N / K);
    if numBlocks < 1
        warning('Not enough samples for K=%d; skipping %s tree.', K, labelName);
        return
    end
    trimmed   = s(1:numBlocks*K);
    blocksMat = reshape(trimmed, K, numBlocks)';                 % [numBlocks x K]
    blkSyms   = cellstr(join(string(blocksMat), "_", 2));

    % Run baseline to get codebook over blocks
    R = baseline_huffman_V2(string(blkSyms));
    codebook = local_extract_codebook(R);
    if isempty(codebook), error('Codebook not found in baseline output.'); end

    % Build trie graph
    [G, nodeLabels, edgeLabels] = local_build_trie_graph(codebook);

    % Plot (dark mode)
    f = figure('Name',sprintf('Huffman Tree — %s (K=%d)', labelName, K), 'Color','k');
    ax = axes('Parent', f, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); %#ok<NASGU>
    p = plot(G, 'Layout','layered', 'Direction','right', 'NodeLabel', nodeLabels); hold on;
    set(gca,'Color','k','XColor','w','YColor','w','Visible','off');
    p.NodeColor      = [0.8 0.8 0.8];
    p.NodeLabelColor = [1   1   1  ];
    p.EdgeColor      = [0.7 0.7 0.7];

    % Edge labels (0/1) at midpoints
    midsX = 0.5*(p.XData(G.Edges.EndNodes(:,1)) + p.XData(G.Edges.EndNodes(:,2)));
    midsY = 0.5*(p.YData(G.Edges.EndNodes(:,1)) + p.YData(G.Edges.EndNodes(:,2)));
    for e = 1:numedges(G)
        text(midsX(e), midsY(e), [' ' edgeLabels{e}], ...
             'Color','w','FontWeight','bold','HorizontalAlignment','left');
    end

    title(sprintf('Huffman Tree — %s (K = %d)', labelName, K), 'Color','w','FontWeight','bold');
end

function codebook = local_extract_codebook(R)
    % Return Nx2 cell array: {symbol, code}
    codebook = {};
    cands = ["codebook","dict","huffman_codebook","huff_codebook","codes"];
    for f = cands
        if isfield(R, f)
            cb = R.(f);
            if iscell(cb) && size(cb,2) >= 2
                sy = string(cb(:,1)); co = string(cb(:,2));
                codebook = [cellstr(sy) cellstr(co)];
                return
            elseif isstruct(cb)
                symF  = ["symbol","sym","token","key","name"];
                codeF = ["code","bits","word"];
                for sf = symF
                    for cf = codeF
                        if all(isfield(cb, [sf cf]))
                            sy = string({cb.(sf)}'); co = string({cb.(cf)}');
                            codebook = [cellstr(sy) cellstr(co)];
                            return
                        end
                    end
                end
            end
        end
    end
    if isfield(R,'symbols') && isfield(R,'codes')
        sy = string(R.symbols(:)); co = string(R.codes(:));
        if numel(sy)==numel(co), codebook = [cellstr(sy) cellstr(co)]; end
    end
end

function [G, nodeLabels, edgeLabels] = local_build_trie_graph(codebook)
    % Build a prefix trie graph from {symbol, code}
    nodeId = containers.Map('KeyType','char','ValueType','int32');
    nodes = {""}; nodeId('') = 1;
    edges_s = []; edges_t = []; edgeLabels = {};
    leafNames = containers.Map('KeyType','char','ValueType','char');

    for i = 1:size(codebook,1)
        sym  = string(codebook{i,1});
        code = char(string(codebook{i,2}));
        cur  = '';
        for k = 1:numel(code)
            bit = code(k);
            next = [cur bit];
            if ~isKey(nodeId, next)
                nodes{end+1} = next; %#ok<AGROW>
                nodeId(next) = numel(nodes);
            end
            s = nodeId(cur); t = nodeId(next);
            edges_s(end+1) = s; %#ok<AGROW>
            edges_t(end+1) = t; %#ok<AGROW>
            edgeLabels{end+1} = bit; %#ok<AGROW>
            cur = next;
        end
        leafNames(cur) = char(sym);
    end

    nodeLabels = nodes;
    for i = 1:numel(nodes)
        key = nodes{i};
        if isKey(leafNames, key)
            nodeLabels{i} = sprintf('%s', leafNames(key));
        else
            if key == ""
                nodeLabels{i} = "root";
            else
                nodeLabels{i} = sprintf('%s', key);
            end
        end
    end
    G = graph(edges_s, edges_t);
end

