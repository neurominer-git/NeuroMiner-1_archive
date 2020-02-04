function CC = cell2matpadnan(C,pad)

if ~exist('pad','var') || isempty(pad), pad=NaN; end 
maxLengthCell=max(cellfun('size',C,2));  %finding the longest vector in the cell array
CC = cell(size(C));
for i=1:length(C)
    CC{i} = [C{i} repmat(pad, 1, maxLengthCell - numel(C{i})) ];   %zeropad the elements in each cell array with a length shorter than the maxlength
end
CC=cell2mat(CC); %A is your matrix