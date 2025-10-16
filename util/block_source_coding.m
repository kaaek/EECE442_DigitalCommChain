function results = block_source_coding(seq, A, K_values, verify_lossless)
% BLOCK_SOURCE_CODING  (e–h) Block source coding with fixed-length vs. Huffman,
% implemented via your baseline_huffman_V2 (no toolbox calls).
%
% Usage:
%   results = block_source_coding(seq, A, [2 3 4], true);

if nargin < 3 || isempty(K_values), K_values = [2 3 4]; end
if nargin < 4, verify_lossless = true; end

seq = seq(:)';                 % row vector
A   = unique(A(:)');           % clean alphabet
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
  blocksMat = reshape(trimmed, K, numBlocks)';                 % [numBlocks x K]

  % Represent each block as a single symbol: e.g., "a_b_c"
  blkSyms = cellstr(join(string(blocksMat), "_", 2));          % 1 x numBlocks (cellstr)

  % Unique block symbols and empirical probabilities p(a^K)
  [uniqBlkSyms, ~, idx] = unique(blkSyms, 'stable');
  counts = accumarray(idx, 1);
  p = counts / sum(counts);

  % ---------- (f)(i) fixed-length baseline ----------
  % We keep the same fixed-length baseline as your original: K * ceil(log2(M))
  % (i.e., symbol-wise fixed-length applied K times per block).
  fixed_bits_per_block  = K * ceil(log2(max(M,1)));
  fixed_bits_per_symbol = fixed_bits_per_block / K;

  % ---------- (f)(ii) Huffman on p(a^K) via baseline_huffman_V2 ----------
  % Call your baseline with the *sequence of block symbols*.
  % We DO NOT pass A_design, because your block fixed-length baseline
  % is defined separately above (not the dict's "designed" alphabet).
  Rblk = baseline_huffman_V2(string(blkSyms));   % uses your build/encode/decode & verifies lossless

  % Average Huffman length per *block* (from your Rblk)
  L_emp_block = Rblk.huffman_avg_bits_per_symbol;           % "symbol" here == one block
  L_emp_sym   = L_emp_block / K;

  % Actual length on THIS data (bitstring actually produced by your encoder)
  total_bits_actual = double(Rblk.total_bits_huffman);      % strlength returns double-compatible
  L_actual_sym = (total_bits_actual / numBlocks) / K;

  % Optional: verify lossless reconstruction (already done inside baseline; re-check flag)
  if verify_lossless
    if ~isfield(Rblk, 'lossless_verified') || ~Rblk.lossless_verified
      error('baseline_huffman_V2 failed lossless verification for K=%d.', K);
    end
    % Rebuild the original trimmed sequence from decoded blocks if desired.
    % Not strictly necessary since baseline already verified, so we skip a second decode.
  end

  % ---------- (g) block entropy H_K and rate H_K/K ----------
  % Use the empirical block distribution p over uniqBlkSyms
  Hk = -sum(p .* log2(p + (p==0)));
  Hk_per_sym = Hk / K;

  gain_emp    = fixed_bits_per_symbol - L_emp_sym;
  gain_actual = fixed_bits_per_symbol - L_actual_sym;

  fprintf('K = %d | blocks = %d | unique blocks = %d\n', K, numBlocks, numel(uniqBlkSyms));
  fprintf('  Fixed length    : %.4f bits/symbol\n', fixed_bits_per_symbol);
  fprintf('  Huffman (emp)   : %.4f bits/symbol   [H_K/K = %.4f]\n', L_emp_sym, Hk_per_sym);
  fprintf('  Huffman (actual): %.4f bits/symbol\n', L_actual_sym);
  fprintf('  Gains vs fixed  : +%.4f (emp), +%.4f (actual)\n\n', gain_emp, gain_actual);

  % Fill results row
  T{r,"K"}                         = K;
  T{r,"NumBlocks"}                 = numBlocks;
  T{r,"Acard"}                     = M;
  T{r,"Fixed_bits_per_sym"}        = fixed_bits_per_symbol;
  T{r,"Huff_bits_per_sym_emp"}     = L_emp_sym;
  T{r,"Huff_bits_per_sym_actual"}  = L_actual_sym;
  T{r,"Hk"}                        = Hk;
  T{r,"Hk_per_sym"}                = Hk_per_sym;
  T{r,"Gain_vs_fixed_emp"}         = gain_emp;
  T{r,"Gain_vs_fixed_actual"}      = gain_actual;
end

results = T;

% ---------- (h) discussion tips ----------
fprintf('--- DISCUSSION HINTS (h) ---\n');
fprintf('* Huffman < fixed-length -> skew/structure in block histogram.\n');
fprintf('* L_emp_sym close to H_K/K -> near-optimal for empirical model.\n');
fprintf('* Larger K helps with short-range dependencies; sparsity limits gains.\n');
fprintf('* Dead-zones / clustered outputs -> bigger gains; uniform histograms -> smaller gains.\n');
end
