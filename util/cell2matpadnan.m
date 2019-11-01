function C = cell2matpadnan(C,pad)

if ~exist('pad','var') || isempty(pad), pad=NaN; end 
maxLengthCell=max(cellfun('size',C,2));  %finding the longest vector in the cell array
for i=1:length(C)
    for j=cellfun('size',C(i),2)+1:maxLengthCell
         C{i}(j)=pad;   %zeropad the elements in each cell array with a length shorter than the maxlength
    end
end
C=cell2mat(C); %A is your matrix