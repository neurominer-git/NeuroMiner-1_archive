function PP = nk_CheckFileExist(P)

PP=[];

if iscell(P), P = char(P);end

for i=1:size(P,1)
    iP = deblank(P(i,:));
    if ~exist(iP,'file')
        if isempty(PP)
            PP = iP;
        else
            PP = char(PP,iP);
        end
    end
    
end