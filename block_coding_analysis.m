% ----------------------------------------------------------------------
% authors: Khalil El Kaaki, Mouhammad Kandakji
%
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% * Dark-mode edits by ChatGPT
% ----------------------------------------------------------------------
function [cmp_tbl, results_uniform, results_lloydmax, decoded_uniform, decoded_lloydmax, is_lossless_uniform, is_lossless_lloydmax] = ...
    block_coding_analysis(seq_uniform, seq_lloydmax, A, K_values, verify_lossless)


fprintf('========================================\n');
fprintf('   3.2 Block Source Coding is Running...\n');
fprintf('========================================\n');

% --------- Input defaults ---------
if nargin < 4 || isempty(K_values), K_values = [1 2 3 4]; end
if nargin < 5, verify_lossless = true; end
if nargin < 2 || isempty(seq_lloydmax), seq_lloydmax = seq_uniform; end

% ======= DEBUG/SAFETY GUARD: verify inputs are correct =======
if islogical(K_values) || (isscalar(K_values) && (K_values==0 || K_values==1))
    error(['block_coding_analysis: Bad K_values detected.\n' ...
           'Call with FIVE args and pass alphabet VALUES.\n' ...
           'Example:\nA_vals = unique(xq_u(:)).''; block_coding_analysis(xq_u, [], A_vals, [1 2 3 4], true);']);
end

% Ensure A is VALUES. If empty or scalar, rebuild from sequences.
if isempty(A) || isscalar(A)
    if isempty(seq_lloydmax)
        A = unique(seq_uniform(:)).';
    else
        A = unique([seq_uniform(:); seq_lloydmax(:)]).';
    end
end
A = unique(A(:)).';

% Normalize K_values shape
K_values = K_values(:).';
if isempty(K_values), K_values = [1 2 3 4]; end

% Echo effective inputs
fprintf('Effective inputs:\n');
fprintf('  |A| (alphabet size): %d\n', numel(A));
fprintf('  K_values: %s\n', mat2str(K_values));
fprintf('  verify_lossless: %d\n', verify_lossless);
fprintf('  N_uniform=%d, N_lloydmax=%d\n\n', numel(seq_uniform), numel(seq_lloydmax));

% Keep originals for lossless checks
seq_syms_u_raw  = seq_uniform(:);
seq_syms_lm_raw = seq_lloydmax(:);

% ===== Run block source coding on BOTH sequences (collect tables + artifacts) =====
results_uniform_all  = cell(numel(K_values),1);
results_lloydmax_all = cell(numel(K_values),1);
artifacts_uniform_all  = cell(numel(K_values),1);
artifacts_lloydmax_all = cell(numel(K_values),1);

for i = 1:numel(K_values)
    K = K_values(i);
    fprintf('Running block_source_coding for K = %d ...\n', K);
    [ru, au]  = block_source_coding(seq_uniform,  A, K, verify_lossless);
    [rlm, alm] = block_source_coding(seq_lloydmax, A, K, verify_lossless);
    results_uniform_all{i}   = ru;
    results_lloydmax_all{i}  = rlm;
    artifacts_uniform_all{i} = au;
    artifacts_lloydmax_all{i}= alm;
end

% ===== Combine per-K TABLE rows into single tables =====
results_uniform  = vertcat(results_uniform_all{:});
results_lloydmax = vertcat(results_lloydmax_all{:});

% ===== Latency metric =====
results_uniform.Latency_symbols  = results_uniform.K;
results_lloydmax.Latency_symbols = results_lloydmax.K;

% ===== Cross-comparison table =====
[Ks_common, iu, il] = intersect(results_uniform.K, results_lloydmax.K, 'stable');
compare_table = {};
for t = 1:numel(Ks_common)
    K = Ks_common(t);
    total_bits_u  = results_uniform.Total_bits(iu(t));
    bps_u         = results_uniform.Bits_per_symbol(iu(t));
    latency_u     = results_uniform.Latency_symbols(iu(t));
    total_bits_lm = results_lloydmax.Total_bits(il(t));
    bps_lm        = results_lloydmax.Bits_per_symbol(il(t));
    latency_lm    = results_lloydmax.Latency_symbols(il(t));
    compare_table(end+1,:) = {"Uniform",   K, total_bits_u,  bps_u,  latency_u}; 
    compare_table(end+1,:) = {"Lloyd–Max", K, total_bits_lm, bps_lm, latency_lm}; 
end
cmp_hdr = {'Quantizer','K','Total_bits','Bits_per_symbol','Latency_symbols'};
cmp_tbl = cell2table(compare_table, 'VariableNames', cmp_hdr);

% ===== Losslessness =====
decoded_uniform  = [];
decoded_lloydmax = [];
is_lossless_uniform  = false;
is_lossless_lloydmax = false;

% find K=1 artifact if present
idxK1_u  = find(results_uniform.K  == 1, 1, 'first');
idxK1_lm = find(results_lloydmax.K == 1, 1, 'first');

if ~isempty(idxK1_u) && ~isempty(artifacts_uniform_all{idxK1_u})
    au = artifacts_uniform_all{idxK1_u};
    if ~isempty(au) && isstruct(au) && ~isempty(au(1).decoded_seq)
        decoded_uniform = au(1).decoded_seq(:);
        is_lossless_uniform = equal_with_tol(decoded_uniform, seq_syms_u_raw(1:numel(decoded_uniform)));
    end
end

if ~isempty(idxK1_lm) && ~isempty(artifacts_lloydmax_all{idxK1_lm})
    al = artifacts_lloydmax_all{idxK1_lm};
    if ~isempty(al) && isstruct(al) && ~isempty(al(1).decoded_seq)
        decoded_lloydmax = al(1).decoded_seq(:);
        is_lossless_lloydmax = equal_with_tol(decoded_lloydmax, seq_syms_lm_raw(1:numel(decoded_lloydmax)));
    end
end

% ===== Display Tables Automatically =====
fprintf('\n==================== UNIFORM RESULTS ====================\n');
disp(results_uniform);
fprintf('\n==================== LLOYD–MAX RESULTS ==================\n');
disp(results_lloydmax);
fprintf('\n==================== UNIFORM vs LLOYD–MAX COMPARISON ====================\n');
disp(cmp_tbl);

% ===== Losslessness summary =====
if verify_lossless
    fprintf('\n==================== LOSSLESSNESS CHECKS (from artifacts, K=1) ====================\n');
    fprintf('Uniform:   %s\n',   ternary(is_lossless_uniform,  'PASS (decoded == original)', ...
                                                        ternary(isempty(decoded_uniform),'N/A (decoded not found)','FAIL (decoded ~= original)')));
    fprintf('Lloyd–Max: %s\n',   ternary(is_lossless_lloydmax, 'PASS (decoded == original)', ...
                                                        ternary(isempty(decoded_lloydmax),'N/A (decoded not found)','FAIL (decoded ~= original)')));
end

% ===== Plot 1: Bits per symbol vs K =====
figure('Name','Bits per symbol vs K','Color','k');
ax1 = axes; set(ax1,'Color','k','XColor','w','YColor','w'); hold on; grid on; ax1.GridColor=[.5 .5 .5];
plot(results_uniform.K, results_uniform.Bits_per_symbol, '-o','Color',[0 .6 1],'LineWidth',1.6,'DisplayName','Uniform');
plot(results_lloydmax.K, results_lloydmax.Bits_per_symbol, '--s','Color',[1 .3 .3],'LineWidth',1.6,'DisplayName','Lloyd–Max');
xlabel('Block size K','Color','w'); ylabel('bits / symbol','Color','w');
title('Throughput (bits per symbol) vs K','Color','w','FontWeight','bold');
legend('TextColor','w','Color',[.1 .1 .1],'Location','best');

% ===== Plot 2: Dependency gains vs K (Uniform) =====
% Pull directly from the table columns
Hk_per_sym       = results_uniform.Hk_per_sym;
L_emp_per_sym    = results_uniform.Huff_bits_per_sym_emp;
dep_gain_emp     = L_emp_per_sym(1) - L_emp_per_sym;   % (ℓ̄1 - ℓ̄K)
dep_gain_th      = Hk_per_sym(1)    - Hk_per_sym;      % (H1 - HK/K)

figure('Name','Dependency Gains','Color','k');
ax2 = axes; set(ax2,'Color','k','XColor','w','YColor','w'); hold on; grid on; ax2.GridColor=[.5 .5 .5];
plot(results_uniform.K, dep_gain_emp, '-o','LineWidth',1.6,'Color',[0 0.6 1],'DisplayName','$\bar{\ell}_{1}-\bar{\ell}_{K}$');
plot(results_uniform.K, dep_gain_th,  '--s','LineWidth',1.6,'Color',[1 0.8 0.2],'DisplayName','$H_{1}-H_{K}/K$');
xlabel('Block size K','Color','w'); ylabel('bits / symbol','Color','w');
title('Gains from Short-Range Dependencies (Uniform)','Color','w','FontWeight','bold');
legend('Interpreter','latex','TextColor','w','Color',[.1 .1 .1],'Location','best');

% ===== Summary printout =====
fprintf('\n==================== DISCUSSION SUMMARY ====================\n');
for i = 1:numel(results_uniform.K)
    fprintf('K=%d | H_K/K=%.4f | l̄_emp=%.4f | Gain(emp)=+%.4f | Gain(th)=+%.4f\n', ...
        results_uniform.K(i), Hk_per_sym(i), L_emp_per_sym(i), ...
        dep_gain_emp(i), dep_gain_th(i));
end

fprintf('\n=== All plots generated automatically. ===\n');
fprintf('========================================\n');

end

% ===================== helpers =====================

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function tf = equal_with_tol(x, y)
    if isempty(x) || isempty(y) || numel(x) ~= numel(y), tf = false; return; end
    try
        xd = double(x(:)); yd = double(y(:));
        tol = 1e-12 * max(1, max(abs([xd; yd])));
        tf = all(abs(xd - yd) <= tol);
    catch
        
        xs = normalize_seq_for_compare(x); ys = normalize_seq_for_compare(y);
        tf = isequal(xs, ys);
    end
end

function y = normalize_seq_for_compare(x)
if isempty(x), y = []; return; end
if isnumeric(x), y = string(compose('%.15g', double(x(:)))); return; end
if iscategorical(x), y = string(x(:)); y = strtrim(y); return; end
if isstring(x), y = strtrim(x(:)); return; end
if ischar(x), y = string(cellstr(x(:))); y = strtrim(y); return; end
if iscellstr(x), y = string(x(:)); y = strtrim(y); return; end
y = string(x(:)); y = strtrim(y);
end
