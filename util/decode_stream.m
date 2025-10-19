% ----------------------------------------------------------------------
% author: Mouhamad Kandakji
% 
% Note on the use of AI:
% * Copilot wrote the help sections for our functions
%       (the big comment blocks following function declarations)
% * ChatGPT only corrected minor logical and syntax errors.
% ----------------------------------------------------------------------

function S = decode_stream(bits, dict)
trie = struct('next', containers.Map({'0','1'},{[],[]}), 'sym', "");
triePool = trie; % root
function id = newNode()
    triePool(end+1) = struct('next', containers.Map({'0','1'},{[],[]}), 'sym', ""); %#ok<AGROW>
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
        out(end+1,1) = triePool(node).sym; %#ok<AGROW>
        node = 1;
    end
end
S = out;
end