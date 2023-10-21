function decoder = treedecoder(code, option)

%treedecoder - Builds a tree decoder for a prefix code.
%
%   decoder = treedecoder(code [, option]) returns the tree decoder of code 
%   'code'.  Assumes the code is prefix.  Useful for passing to the function 
%   'decode' so that the tree decoder does not has to be generated at each 
%   decoding.
%
%   Code 'code' must be a cell of char (row) arrays.  Each element of this array
%   is the binary codeword.
%
%   Hashmarks are shown during the construction of the lookup table if the option 
%   'verbose' is passed.
%
%   See also decode, lookupdecode, treedecode, lookupdecoder.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~iscell(code) | ~isrow(code)
    error('code must be a row cell vector');
end

if nargin < 2
    verbose = false;
else
    if ~isa(option, 'char')
        error('option must be a string');
    end
    verbose = strcmp(option, 'verbose');
end

emptynode.branches = [0 0];
emptynode.index = 0;

if(length(code) == 0)
    tree = emptynode(ones(1, 0));
    return
end

tree(1) = emptynode;

n_symbols = 0;

hash = ceil(length(code) / 79);

for i = 1 : length(code)
    if verbose
        if mod(i - 1, hash) == 0
            fprintf('#');
        end
    end
    
    if(code{i} == '!')
        continue;
    end

    n_symbols = n_symbols + 1;
    
    codeword = logical(code{i}); 
    
    % Find node corresponding to codeword:
    node = 1;
    for j = 1 : length(codeword)
        % Follow branch according to bit, build node if necessary.
        if(tree(node).index ~= 0)
            error(['non prefix code: ' code{tree(node).index} ' and ' code{i}]);
        end
        if(tree(node).branches(codeword(j) + 1) == 0)
            new = length(tree) + 1;
            tree(new) = emptynode;
            tree(node).branches(codeword(j) + 1) = new;
        end
        node = tree(node).branches(codeword(j) + 1);
    end

    % Store index of symbol in tree:
    if(~isemptynode(tree(node)))
        error(['non prefix code: ' code{i} 'is prefix of others']);
    end
    tree(node).index = i;
end

if verbose
    fprintf('\n');
end

redundant = false;
i = 1;
while(i <= length(tree) & (isfullnode(tree(i)) | isemptynode(tree(i))))
    i = i + 1;
end
if(i <= length(tree))
    warning('code has redundant codewords');
end

decoder.tree = tree;
decoder.n_symbols = n_symbols;
decoder.type = 'tree';
decoder.decode = @treedecode;

function r = isfullnode(n)

r = all(n.branches ~= [0 0]);

function r = isemptynode(n)

r = all(n.branches == [0 0]);
