function R = baseline_huffman(x)
% BASELINE_HUFFMAN  Baseline Huffman Coding (symbol-wise) + Tree Plot
% x : vector of discrete symbols (numeric, char, or strings)
% R : struct with fields shown at the bottom

% ---------- normalize symbols ----------
if isstring(x) || ischar(x); x = cellstr(string(x(:))); end
if iscell(x)
    S = string(x(:));
else
    S = string(x(:));
end
N = numel(S);

% ---------- (a) Unknown distribution pass: fixed-length baseline ----------
A = unique(S,'stable');                 % alphabet (stable order)
M = numel(A);
L_fixed_per_sym = ceil(log2(max(M,1))); % handle M=1 safely
bits_fixed = N * L_fixed_per_sym;

% ---------- (b) Empirical modeling & entropy ----------
[~,~,idx] = unique(S,'stable');   % FIXED: use 'stable' to align A & p
cnt = accumarray(idx,1,[M,1]);
p = cnt / N;
H = -sum(p .* log2(p + (p==0)));        % avoid log2(0)

% ---------- (c) Optimal instantaneous (Huffman) code ----------
[dict, root] = build_huffman_dict(A, p);   % returns dictionary & root

% Map symbol -> probability (for dict alignment)
sym2p = containers.Map(cellstr(A), num2cell(p));

% Compute average code length correctly
L_avg = 0;
codeLen   = zeros(size(dict,1),1);
p_aligned = zeros(size(dict,1),1);
for i = 1:size(dict,1)
    ai = char(dict{i,1});     % symbol
    ci = char(dict{i,2});     % code
    li = strlength(ci);       % code length
    pi = sym2p(ai);           % correct p(a)
    codeLen(i)   = li;
    p_aligned(i) = pi;
    L_avg        = L_avg + double(li) * pi;
end

% Encode
codeMap = containers.Map(cellstr(dict(:,1)), cellstr(dict(:,2)));
encoded = encode_stream(S, codeMap);

% Decode and verify
decoded = decode_stream(encoded, dict);
ok = all(S == string(decoded));

% ---------- Print each unique symbol exactly once ----------
fprintf('\n--- Huffman Code Dictionary (unique symbols) ---\n');
fprintf(' symbol           p(a)      code     len\n');
for i = 1:size(dict,1)
    sym = char(dict{i,1});
    pj  = p_aligned(i);
    cj  = char(dict{i,2});
    fprintf(' %-15s  %7.4f    %-8s %d\n', sym, pj, cj, strlength(cj));
end

% ---------- (d) Analysis ----------
R = struct();
R.N = N;
R.alphabet = A;
R.counts = cnt;
R.p = p;
R.entropy_bits_per_symbol = H;
R.fixed_bits_per_symbol = L_fixed_per_sym;
R.huffman_avg_bits_per_symbol = L_avg;
R.total_bits_fixed = bits_fixed;
R.total_bits_huffman = numel(encoded);
R.compression_gain_vs_fixed = 1 - R.total_bits_huffman / bits_fixed;
R.lossless_verified = ok;
R.dict = dict;
R.encoded_bitstring = encoded;          % as char('0'/'1')
R.dict_table = table(string(dict(:,1)), p_aligned(:), string(dict(:,2)), double(codeLen(:)), ...
    'VariableNames', {'Symbol','Probability','Code','Len'});
R.tree_root = root;

% Sanity checks
assert(L_avg + 1e-12 >= H, 'Avg code length should be >= entropy');
assert(L_avg <= L_fixed_per_sym + 1e-12, 'Avg code length should be <= fixed length');

% ---------- Draw the Huffman tree ----------
figure('Name','Huffman Tree','Color','w');
plot_huffman_tree(root);
title('Huffman Code Tree');

% Pretty print summary
fprintf('\n--- Baseline Huffman Report ---\n');
fprintf('N=%d, |A|=%d\n', N, M);
fprintf('Entropy H(A)          = %.4f bits/sym\n', H);
fprintf('Fixed length          = %d bits/sym (total %d)\n', L_fixed_per_sym, bits_fixed);
fprintf('Huffman avg length    = %.4f bits/sym (total %d)\n', L_avg, numel(encoded));
fprintf('Compression vs fixed  = %.2f%%\n', 100*R.compression_gain_vs_fixed);
fprintf('Lossless verified     = %d\n', ok);
end

% ===== Helpers =====
function [dict, root] = build_huffman_dict(A, p)
nodes = struct('sym',[], 'p',[], 'left',[], 'right',[]);
K = numel(A);
T = repmat(nodes, K, 1);
for i = 1:K
    T(i).sym = string(A(i));
    T(i).p   = p(i);
    T(i).left = [];
    T(i).right = [];
end
forest = num2cell(T(:));

while numel(forest) > 1
    [~,order] = sort(cellfun(@(n) n.p, forest)); % ascending prob
    forest = forest(order);
    a = forest{1}; 
    b = forest{2};
    parent = struct('sym',"", 'p', a.p + b.p, 'left', a, 'right', b);
    forest = [{parent}; forest(3:end)];
end
root = forest{1};

pairs = traverse(root, "");
dict = cell(numel(pairs), 2);
for i = 1:numel(pairs)
    dict{i,1} = pairs(i).sym;
    dict{i,2} = pairs(i).code;
end
end

function out = traverse(node, prefix)
if isempty(node.left) && isempty(node.right)
    if prefix == ""      % edge case: single-symbol alphabet
        prefix = "0";
    end
    out = struct('sym', node.sym, 'code', char(prefix));
    return
end
out = struct('sym', {}, 'code', {});
if ~isempty(node.left)
    out = [out, traverse(node.left, prefix + "0")]; 
end
if ~isempty(node.right)
    out = [out, traverse(node.right, prefix + "1")]; 
end
end

function bits = encode_stream(S, codeMap)
buf = strings(numel(S),1);
for i = 1:numel(S)
    buf(i) = string(codeMap(char(S(i))));
end
bits = char(strjoin(buf, ""));
end

function S = decode_stream(bits, dict)
trie = struct('next', containers.Map({'0','1'},{[],[]}), 'sym', "");
triePool = trie; % root
function id = newNode()
    triePool(end+1) = struct('next', containers.Map({'0','1'},{[],[]}), 'sym', ""); 
    id = numel(triePool);
end
for i = 1:size(dict,1)
    sym = string(dict{i,1}); code = char(dict{i,2});
    node = 1;
    for c = code
        nxt = triePool(node).next(c);
        if isempty(nxt)
            nxt = newNode();
            triePool(node).next(c) = nxt;
        end
        node = nxt;
    end
    triePool(node).sym = sym;
end
out = strings(0,1);
node = 1;
for c = bits
    node = triePool(node).next(c);
    if triePool(node).sym ~= ""
        out(end+1,1) = triePool(node).sym; 
        node = 1;
    end
end
S = out;
end

function plot_huffman_tree(root)
nodes = struct('sym',{},'p',{},'left',{},'right',{});
edges = [];
function id = addNode(node)
    id = numel(nodes)+1;
    nodes(id).sym   = node.sym;
    nodes(id).p     = node.p;
    nodes(id).left  = [];
    nodes(id).right = [];
    if ~isempty(node.left)
        L = addNode(node.left);
        nodes(id).left = L;
        edges(end+1,:) = [id, L, 0]; 
    end
    if ~isempty(node.right)
        R = addNode(node.right);
        nodes(id).right = R;
        edges(end+1,:) = [id, R, 1]; 
    end
end
rootId = addNode(root);

n = numel(nodes);
pos = zeros(n,2);
leafX = 0;
function assignPos(id, depth)
    isLeaf = isempty(nodes(id).left) && isempty(nodes(id).right);
    if isLeaf
        leafX = leafX + 1;
        pos(id,:) = [leafX, depth];
    else
        if ~isempty(nodes(id).left),  assignPos(nodes(id).left,  depth+1); end
        if ~isempty(nodes(id).right), assignPos(nodes(id).right, depth+1); end
        xs = [];
        if ~isempty(nodes(id).left),  xs(end+1) = pos(nodes(id).left,1);  end 
        if ~isempty(nodes(id).right), xs(end+1) = pos(nodes(id).right,1); end 
        pos(id,:) = [mean(xs), depth];
    end
end
assignPos(rootId, 0);
set(gca,'YDir','reverse');

hold on; axis off;
for k = 1:size(edges,1)
    p = edges(k,1); c = edges(k,2); b = edges(k,3);
    plot([pos(p,1) pos(c,1)], [pos(p,2) pos(c,2)], '-');
    mx = (pos(p,1)+pos(c,1))/2; my = (pos(p,2)+pos(c,2))/2;
    text(mx, my, num2str(b), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'Interpreter','tex');
end

for i = 1:n
    isLeaf = isempty(nodes(i).left) && isempty(nodes(i).right);
    plot(pos(i,1), pos(i,2), 'o', 'MarkerFaceColor', [0.92 0.92 0.92], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    if isLeaf && nodes(i).sym ~= ""
        label = sprintf('%s\\newline p=%.3f', char(nodes(i).sym), nodes(i).p);
    else
        label = sprintf('p=%.3f', nodes(i).p);
    end
    text(pos(i,1)+0.02, pos(i,2)-0.05, label, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Interpreter','tex');
end
axis equal
xlim([min(pos(:,1))-0.5, max(pos(:,1))+0.5])
ylim([min(pos(:,2))-0.5, max(pos(:,2))+0.5])
hold off
end
