function [results, artifacts] = block_source_coding(seq, A, K_values, verify_lossless)
% BLOCK_SOURCE_CODING with built-in decode display + bit throughput, robust artifacts

if nargin < 3 || isempty(K_values), K_values = [2 3 4]; end
if nargin < 4, verify_lossless = true; end

seq = seq(:)';                 
A   = unique(A(:)');           
M   = numel(A);
N   = numel(seq);

fprintf('--- BLOCK SOURCE CODING with Decode + Bit Throughput ---\n');
fprintf('Alphabet size |A| = %d, N = %d\n\n', M, N);

% Build table
results = table('Size',[0 12], ...
  'VariableTypes', repmat("double",1,12), ...
  'VariableNames', ["K","NumBlocks","Acard","Fixed_bits_per_sym","Huff_bits_per_sym_emp", ...
    "Huff_bits_per_sym_actual","Hk","Hk_per_sym","Gain_vs_fixed_emp","Gain_vs_fixed_actual", ...
    "Total_bits","Bits_per_symbol"]);

artifacts = struct( ...
  'K', [], 'numBlocks', [], 'uniqBlkSyms', {{}}, 'p', [], ...
  'trimmed', [], 'blocksMat', [], 'blkSyms', {{}}, ...
  'decoded_blks', {{}}, 'decoded_seq', '', ...
  'encoded_stream', '', 'total_bits', NaN, 'bits_per_symbol', NaN, ...
  'Rblk', struct(), ...
  'fixed_bits_per_block', [], 'fixed_bits_per_symbol', [], ...
  'L_emp_block', [], 'L_emp_sym', [], 'total_bits_actual', [], 'L_actual_sym', [], ...
  'Hk', [], 'Hk_per_sym', [] );
artifacts(1) = [];  

for K = K_values
  % ---------- (e) form non-overlapping K-blocks ----------
  numBlocks = floor(N / K);
  if numBlocks < 1
    warning('K=%d too large for N=%d. Skipping.', K, N);
    continue
  end
  trimmed   = seq(1:numBlocks*K);
  blocksMat = reshape(trimmed, K, numBlocks)';
  blkSyms   = cellstr(join(string(blocksMat), "_", 2));

  [uniqBlkSyms, ~, idx] = unique(blkSyms, 'stable');
  counts = accumarray(idx, 1);
  p = counts / sum(counts);

  % ---------- (f)(i) fixed-length baseline ----------
  fixed_bits_per_block  = K * ceil(log2(max(M,1)));
  fixed_bits_per_symbol = fixed_bits_per_block / K;

  % ---------- (f)(ii) Huffman encode/decode ----------
  Rblk = baseline_huffman_V2(string(blkSyms));   

  L_emp_block = fetchfield_safe(Rblk, ["huffman_avg_bits_per_symbol","avg_bits_per_symbol","Lavg_block"], NaN);
  L_emp_sym   = L_emp_block / K;

  total_bits_actual = double(fetchfield_safe(Rblk, ["total_bits_huffman","total_bits","bitcount","n_bits"], NaN));
  L_actual_sym = (total_bits_actual / numBlocks) / K;

  % ---------- Bit count + throughput extraction  ----------
  encoded_stream = fetchfield_safe(Rblk, ["encoded_stream","bitstream","encoded_bits","code"], []);
  if iscell(encoded_stream), encoded_stream = [encoded_stream{:}]; end
  if isstring(encoded_stream), encoded_stream = char(encoded_stream); end
  if isnumeric(encoded_stream), encoded_stream = char(string(encoded_stream)); end
  if isempty(encoded_stream) && ~isnan(total_bits_actual)
      % Fall back to the baseline-reported total
      total_bits = total_bits_actual;
  else
      total_bits = numel(encoded_stream);
  end
  bits_per_symbol = total_bits / numel(blkSyms);
  fprintf('K=%d  Bit count: %d bits, Throughput: %.4f bits/symbol\n', K, total_bits, bits_per_symbol);

  % ---------- Decode + verify ----------
  [decoded_blks, decoded_seq] = local_get_decoded_from_baseline(Rblk, K);
  if verify_lossless
    ok = false;
    if ~isempty(decoded_blks)
      ok = (numel(decoded_blks) == numel(blkSyms)) && all(string(decoded_blks) == string(blkSyms));
    end
    if ~ok && ~isempty(decoded_seq)
      ok = isequal(char(decoded_seq(:).'), char(trimmed(:).'));
    end
    if ~ok
      fprintf('Lossless check failed for K=%d: decoded output does not match input.', K);
    else
      fprintf('  Lossless (baseline decode): OK ✅\n');
    end
  end

  % ---------- Entropy ----------
  Hk = -sum(p .* log2(p + (p==0)));
  Hk_per_sym = Hk / K;

  gain_emp    = fixed_bits_per_symbol - L_emp_sym;
  gain_actual = fixed_bits_per_symbol - L_actual_sym;

  % ---------- Print brief summary ----------
  fprintf('  Fixed=%.4f | Huff(emp)=%.4f | Huff(actual)=%.4f | H_K/K=%.4f | Gain(emp)=%.4f | Gain(actual)=%.4f\n\n', ...
    fixed_bits_per_symbol, L_emp_sym, L_actual_sym, Hk_per_sym, gain_emp, gain_actual);

  % ---------- Append table row (safe) ----------
  results(end+1,:) = {K, numBlocks, M, fixed_bits_per_symbol, L_emp_sym, L_actual_sym, ...
                      Hk, Hk_per_sym, gain_emp, gain_actual, total_bits, bits_per_symbol};

  % ---------- Build artifact with the SAME fixed schema every time ----------
  S = struct( ...
    'K', K, ...
    'numBlocks', numBlocks, ...
    'uniqBlkSyms', {uniqBlkSyms}, ...
    'p', p, ...
    'trimmed', trimmed, ...
    'blocksMat', blocksMat, ...
    'blkSyms', {blkSyms}, ...
    'decoded_blks', {decoded_blks}, ...
    'decoded_seq', char(decoded_seq), ...
    'encoded_stream', char(encoded_stream), ...
    'total_bits', total_bits, ...
    'bits_per_symbol', bits_per_symbol, ...
    'Rblk', Rblk, ...
    'fixed_bits_per_block', fixed_bits_per_block, ...
    'fixed_bits_per_symbol', fixed_bits_per_symbol, ...
    'L_emp_block', L_emp_block, ...
    'L_emp_sym', L_emp_sym, ...
    'total_bits_actual', total_bits_actual, ...
    'L_actual_sym', L_actual_sym, ...
    'Hk', Hk, ...
    'Hk_per_sym', Hk_per_sym );
  artifacts(end+1) = S;     
end

fprintf('--- Done ---\n');
end


% ================= helpers =================
function val = fetchfield_safe(S, names, defaultVal)
val = defaultVal;
for k = 1:numel(names)
  f = names(k);
  if isfield(S, f)
    val = S.(f);
    return
  end
end
end

function [decoded_blks, decoded_seq] = local_get_decoded_from_baseline(Rblk, K)
decoded_blks = {};
decoded_seq  = '';

blk_candidates = ["decoded_blocks","decoded_syms","decoded_symbols","decoded","decoded_block_symbols"];
seq_candidates = ["decoded_sequence","decoded_seq","decoded_chars","decoded_string","reconstructed_sequence"];

for f = blk_candidates
  if isfield(Rblk, f)
    tmp = Rblk.(f);
    if iscell(tmp), decoded_blks = tmp;
    elseif isstring(tmp), decoded_blks = cellstr(tmp);
    end
    if ~isempty(decoded_blks), break; end
  end
end

for f = seq_candidates
  if isfield(Rblk, f)
    tmp = Rblk.(f);
    if isstring(tmp) || ischar(tmp), decoded_seq = char(tmp);
    elseif isnumeric(tmp), decoded_seq = char(tmp);
    end
    if ~isempty(decoded_seq), break; end
  end
end

if isempty(decoded_seq) && ~isempty(decoded_blks)
  pieces = regexp(string(decoded_blks), "_", "split");
  flat = [pieces{:}];
  decoded_seq = char(join([flat{:}], ''));
  decoded_seq = decoded_seq(:).';
end
end
