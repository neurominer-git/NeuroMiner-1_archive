function P = nk_ConvProbabilities(C, ngroups)

m = size(C,1);
P = zeros(m,ngroups);

if iscell(C)

    for i=1:m
        for j=1:ngroups
            P(i,j) = sum(C{i}==j)/numel(C{i});
        end
    end

else
    
    for j=1:ngroups
        P(:,j) = sum(C==j,2) / size( C,2 );
    end
    
end