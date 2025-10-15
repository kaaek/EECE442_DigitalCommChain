function results = block_source_coding(seq, A, K_values, verify_lossless)
% BLOCK_SOURCE_CODING  Parts (e–h) — block coding & analysis (separate file).
%
% Usage:
%   results = block_source_coding(seq, A, [2 3 4], true);
%
% Inputs:
%   seq            : vector of quantized symbols (row or col)
%   A              : alphabet (vector of all possible quantizer output symbols)
%   K_values       : vector of block lengths to test (e.g., [2 3 4])
%   verify_lossless: logical; if true, encodes + decodes to check exact recovery
%
% Outputs:
%   results: table with columns:
%       K, NumBlocks, |A|, Fixed_bits_per_sym, Huffman_bits_per_sym_emp,
%       Huffman_bits_per_sym_actual, Hk, Hk_per_sym, Gain_vs_fixed_emp,
%       Gain_vs_fixed_actual
%
% Notes:
% - Non-overlapping blocks; any tail < K is dropped.
% - Fixed length uses k * ceil(log2 |A|) bits per block.
% - Huffman dictionary is built on *observed* blocks with empirical probs.
% - *_emp uses expected length (sum p * len(codeword)); *_actual measures the
%   length of the actual encoded bitstream of your data.
%
% Requires Communications Toolbox (huffmandict, huffmanenco, huffmandeco).

% --------- defaults & checks ---------
if nargin < 3 || isempty(K_values), K_values = [2 3 4]; end
if nargin < 4, verify_lossless = true; end

seq = seq(:)';                     % ensure row
A   = unique(A(:)');               % sanitize alphabet
M   = numel(A);
N   = numel(seq);

if N < min(K_values)
    error('Sequence too short for the requested block sizes.');
end

% map symbols to a compact label space if needed (robust to arbitrary A)
[~, loc] = ismember(seq, A);
if any(loc == 0)
    error('Some symbols in seq are not in A.');
end

% --------- results store ---------
Rows = numel(K_values);
T = table('Size',[Rows 10], ...
    'VariableTypes', ["double","double","double","double","double","double","double","double","double","double"], ...
    'VariableNames', ["K","NumBlocks","Acard","Fixed_bits_per_sym","Huff_bits_per_sym_emp", ...
                      "Huff_bits_per_sym_actual","Hk","Hk_per_sym","Gain_vs_fixed_emp","Gain_vs_fixed_actual"]);

r = 0;

fprintf('--- BLOCK SOURCE CODING (e–h) ---\n');
fprintf('Alphabet size |A| = %d\n\n', M);

for K = K_values
    r = r + 1;

    %% (e) Form length-K non-overlapping blocks
    numBlocks = floor(N / K);
    trimmed   = seq(1 : numBlocks*K);
    blocksMat = reshape(trimmed, K, numBlocks)';               % [numBlocks x K]

    % Convert each block to a symbol (string) so Huffman can work on block-tokens
    blkSyms = join(string(blocksMat), "_", 2);                  % string column
    blkSyms = cellstr(blkSyms);                                 % cellstr column for dict/enco

    % Empirical distribution p(a^k)
    [uniqBlkSyms, ~, idx] = unique(blkSyms, 'stable');
    counts = accumarray(idx, 1);
    p = counts / sum(counts);

    %% (f) (i) Fixed-length coding for blocks
    fixed_bits_per_block  = K * ceil(log2(M));
    fixed_bits_per_symbol = fixed_bits_per_block / K;

    %% (f) (ii) Huffman on p(a^k)
    % Build dictionary on observed blocks only
    dict = huffmandict(uniqBlkSyms, p);

    % Expected (empirical) bits per block from codeword lengths:
    codeLens = cellfun(@numel, dict(:,2));
    % Align lengths to probabilities in same order as uniqBlkSyms
    [~, pos] = ismember(dict(:,1), uniqBlkSyms);
    L_emp_block = sum(p(pos) .* codeLens);
    L_emp_sym   = L_emp_block / K;

    % Actual bits for THIS stream
    bitstream = huffmanenco(blkSyms, dict);              % vector of 0/1 (double)
    L_actual_blockstream = numel(bitstream);
    L_actual_sym = (L_actual_blockstream / numBlocks) / K;

    %% Optional: lossless verification
    if verify_lossless
        decodedBlkSyms = huffmandeco(bitstream, dict);
        if ~isequal(decodedBlkSyms(:), blkSyms(:))
            error('Huffman decode != encode (block symbol mismatch).');
        end
        % Rebuild sequence and compare
        recBlocks = split_block_symbols(decodedBlkSyms, K);  % [numBlocks x K] double
        recSeq = reshape(recBlocks', 1, []);                 % row vector
        if ~isequal(recSeq, trimmed)
            error('Huffman decode != original trimmed sequence.');
        end
    end

    %% (g) Entropy H_k and H_k/k
    Hk = -sum(p .* log2(p));
    Hk_per_sym = Hk / K;

    %% (g) Report & (h) quick commentary
    gain_emp    = fixed_bits_per_symbol - L_emp_sym;
    gain_actual = fixed_bits_per_symbol - L_actual_sym;

    fprintf('K = %d | blocks = %d | unique blocks = %d\n', K, numBlocks, numel(uniqBlkSyms));
    fprintf('  Fixed length    : %.4f bits/symbol\n', fixed_bits_per_symbol);
    fprintf('  Huffman (emp)   : %.4f bits/symbol   [H_k/k = %.4f]\n', L_emp_sym, Hk_per_sym);
    fprintf('  Huffman (actual): %.4f bits/symbol\n', L_actual_sym);
    fprintf('  Gains vs fixed  : +%.4f (emp), +%.4f (actual)\n\n', gain_emp, gain_actual);

    % Fill results table
    T{r,"K"}                        = K;
    T{r,"NumBlocks"}               = numBlocks;
    T{r,"Acard"}                   = M;
    T{r,"Fixed_bits_per_sym"}      = fixed_bits_per_symbol;
    T{r,"Huff_bits_per_sym_emp"}   = L_emp_sym;
    T{r,"Huff_bits_per_sym_actual"}= L_actual_sym;
    T{r,"Hk"}                      = Hk;
    T{r,"Hk_per_sym"}              = Hk_per_sym;
    T{r,"Gain_vs_fixed_emp"}       = gain_emp;
    T{r,"Gain_vs_fixed_actual"}    = gain_actual;
end

results = T;

% ---- (h) Discussion helper (printed hints) ----
fprintf('--- DISCUSSION HINTS (h) ---\n');
fprintf('* If Huffman bits/sym < fixed length, blocks exploit skew/repetition.\n');
fprintf('* If L_emp_sym ≈ H_k/k, your code is near optimal for the observed model.\n');
fprintf('* Larger K can help more when short-range dependencies exist, but sparsity\n');
fprintf('  (few repeats of many block types) limits gains and can even hurt robustness.\n');
fprintf('* Quantizers with dead-zones or clustered outputs often yield skewed block\n');
fprintf('  histograms -> bigger gains. Memoryless/uniform histograms show smaller gains.\n');

end % main function


% ---------- helpers ----------
function blocks = split_block_symbols(blkSyms, K)
% Convert cellstr of 'a_b_...' back into numeric matrix [numBlocks x K]
numBlocks = numel(blkSyms);
blocks = zeros(numBlocks, K);
for i = 1:numBlocks
    parts = split(blkSyms{i}, '_');
    blocks(i,:) = str2double(parts(:)).';
end
end
