function results = block_source_coding(seq, A, K_values, verify_lossless)
% BLOCK_SOURCE_CODING  (e–h) Block source coding with fixed-length vs. Huffman.
% Usage:
%   results = block_source_coding(seq, A, [2 3 4], true);

if nargin < 3 || isempty(K_values), K_values = [2 3 4]; end
if nargin < 4, verify_lossless = true; end

seq = seq(:)';                % row
A   = unique(A(:)');          % clean
M   = numel(A);
N   = numel(seq);

fprintf('--- BLOCK SOURCE CODING (e–h) ---\n');
fprintf('Alphabet size |A| = %d, N = %d\n\n', M, N);

T = table('Size',[numel(K_values) 10], ...
  'VariableTypes', ["double","double","double","double","double","double","double","double","double","double"], ...
  'VariableNames', ["K","NumBlocks","Acard","Fixed_bits_per_sym","Huff_bits_per_sym_emp", ...
    "Huff_bits_per_sym_actual","Hk","Hk_per_sym","Gain_vs_fixed_emp","Gain_vs_fixed_actual"]);

r = 0;

for K = K_values
  r = r + 1;

  % ---------- (e) form non-overlapping K-blocks ----------
  numBlocks = floor(N / K);
  if numBlocks < 1
    warning('K=%d too large for N=%d. Skipping.', K, N);
    continue
  end
  trimmed   = seq(1:numBlocks*K);
  blocksMat = reshape(trimmed, K, numBlocks)';                % [numBlocks x K]

  % Convert each block to a single "symbol" string like "a_b_c"
  blkSyms = cellstr(join(string(blocksMat), "_", 2));         % cellstr, length=numBlocks

  % Unique block symbols and empirical probabilities p(a^k)
  [uniqBlkSyms, ~, idx] = unique(blkSyms, 'stable');
  counts = accumarray(idx, 1);
  p = counts / sum(counts);

  % ---------- (f)(i) fixed-length coding ----------
  fixed_bits_per_block  = K * ceil(log2(M));
  fixed_bits_per_symbol = fixed_bits_per_block / K;

  % ---------- (f)(ii) Huffman on p(a^k) ----------
  % Build dictionary on the observed blocks only
  dict = huffmandict(uniqBlkSyms, p);

  % Expected (empirical) length from codeword lengths
  codeLens = cellfun(@numel, dict(:,2));                       % lengths in bits
  % Map lengths back to uniqBlkSyms order
  [~, pos] = ismember(uniqBlkSyms, dict(:,1));
  L_emp_block = sum(p .* codeLens(pos));
  L_emp_sym   = L_emp_block / K;

  % Actual bit length on THIS data
  bitstream = huffmanenco(blkSyms, dict);                      % vector 0/1
  L_actual_sym = (numel(bitstream) / numBlocks) / K;

  % Optional: verify lossless reconstruction
  if verify_lossless
    decBlkSyms = huffmandeco(bitstream, dict);
    if ~isequal(decBlkSyms(:), blkSyms(:))
      error('Decode != encode for K=%d (block symbol mismatch).', K);
    end
    % Rebuild sequence and compare to trimmed
    recBlocks = zeros(numBlocks, K);
    for i = 1:numBlocks
      parts = split(decBlkSyms{i}, '_');
      recBlocks(i,:) = str2double(parts(:)).';
    end
    recSeq = reshape(recBlocks', 1, []);
    if ~isequal(recSeq, trimmed)
      error('Recovered sequence mismatch for K=%d.', K);
    end
  end

  % ---------- (g) block entropy H_k and rate H_k/k ----------
  Hk = -sum(p .* log2(p));
  Hk_per_sym = Hk / K;

  gain_emp    = fixed_bits_per_symbol - L_emp_sym;
  gain_actual = fixed_bits_per_symbol - L_actual_sym;

  fprintf('K = %d | blocks = %d | unique blocks = %d\n', K, numBlocks, numel(uniqBlkSyms));
  fprintf('  Fixed length    : %.4f bits/symbol\n', fixed_bits_per_symbol);
  fprintf('  Huffman (emp)   : %.4f bits/symbol   [H_k/k = %.4f]\n', L_emp_sym, Hk_per_sym);
  fprintf('  Huffman (actual): %.4f bits/symbol\n', L_actual_sym);
  fprintf('  Gains vs fixed  : +%.4f (emp), +%.4f (actual)\n\n', gain_emp, gain_actual);

  % Fill results row
  T{r,"K"}                         = K;
  T{r,"NumBlocks"}                = numBlocks;
  T{r,"Acard"}                    = M;
  T{r,"Fixed_bits_per_sym"}       = fixed_bits_per_symbol;
  T{r,"Huff_bits_per_sym_emp"}    = L_emp_sym;
  T{r,"Huff_bits_per_sym_actual"} = L_actual_sym;
  T{r,"Hk"}                       = Hk;
  T{r,"Hk_per_sym"}               = Hk_per_sym;
  T{r,"Gain_vs_fixed_emp"}        = gain_emp;
  T{r,"Gain_vs_fixed_actual"}     = gain_actual;
end

results = T;

% ---------- (h) discussion tips ----------
fprintf('--- DISCUSSION HINTS (h) ---\n');
fprintf('* Huffman < fixed-length -> skew/structure in block histogram.\n');
fprintf('* L_emp_sym close to H_k/k -> near-optimal for empirical model.\n');
fprintf('* Larger K helps with short-range dependencies; sparsity limits gains.\n');
fprintf('* Dead-zones / clustered outputs -> bigger gains; uniform histograms -> smaller gains.\n');
end
