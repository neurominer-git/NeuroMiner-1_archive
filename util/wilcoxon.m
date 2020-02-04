function STATS = wilcoxon(x1, x2, alpha, xnames, ynames, filename)
% Based on the WILCOXON code of Giuseppe Cardillo

dff=sort(x2-x1); %difference between x1 and x2
dff(dff==0)=[]; %eliminate null variations
n=length(dff); %number of ranks
if length(x1)~=n %tell me if there are null variations
    fprintf('There are %d null variations that will be deleted\n',length(x1)-n)
end
if isempty(dff) %if all variations are null variations exit function
    error('There are not variations. Wilcoxon test can''t be performed')
end

%Ranks of absolute value of samples differences with sign
[r,t]   = nm_tiedrank(abs(dff)); %ranks and ties
W       = sum(r.*sign(dff)); %Wilcoxon statics (sum of ranks with sign)
pem     = median(dff); %point estimation of median of differences
m       = ceil(n/2); %location of the median

if mod(n,2) == 0 %If the length of the series is even
    tmp = [dff(1:m) pem dff(m:end)]; %add the median in the middle
    dff = tmp; clear tmp 
    m   = m+1;
end

%find how many differences far from the median we have to choose
C       = cumsum(binopdf(0:1:n,n,0.5));
T       = find( C <= alpha/2,1,'last')-1;
cintpem = dff([m-T m+T]); %construct the interval
clear C T m
[I,J]   = ndgrid(dff,dff); d=triu(I+J)./2; %Walsh averages triangular matrix
ld      = sort(d(d~=0)); %linearization of Walsh averages matrix
clear I J 
HLe     = median(ld); %Hodges-Lehmann estimator
if n>15
    A=n*(n+1)/4; B=realsqrt(n*(n+1)*(2*n+1)/24); 
    Za=-realsqrt(2).*erfcinv(2.*(1-alpha/2));
    T=fix(A-Za.*B);
else
    TC=[0 0 0 0 0 0 2 3 5 8 10 13 17 21 25];
    T=TC(n); clear TC
end
cintHLe = ld([T+1 end-T])';
clear dff d ld T%clear unnecessary variable

%If the number of elements N<15 calculate the exact distribution of the
%signed ranks (the number of combinations is 2^N); else use the normal
%distribution approximation.
if n<=15
    ap=ff2n(n); %the all possibilities based on two-level full-factorial design.
    ap(ap~=1)=-1; %change 0 with -1
    k=1:1:n; 
    J=ap*k'; %all possible sums of ranks for k elements
    %to compute the p-value see how many values are more extreme of the observed
    %W and then divide for the total number of combinations
    p=length(J(abs(J)>=abs(W)))/length(J); %p-value
    STATS.method='Exact distribution';
    STATS.W=W;
    STA
else
    sW = sqrt((2*n^3+3*n^2+n-t)/6); %standard deviation
    zW = (abs(W)-0.5)/sW; %z-value with correction for continuity
    p=1-normcdf(zW); %p-value
    STATS.method='Normal approximation';
    STATS.W=W;
    STATS.mean=0;
    STATS.std=sW;
    STATS.z=zW;
    STATS.p=p;
end

if exist('filename','var') && ~isempty(filename)
   sh = {'pp','stats'};
   tbl_perf = array2table([x1 x2],'VariableNames',ynames', 'RowNames', xnames);
   tbl_stats = struct2table(STATS);
   writetable(tbl_perf,filename,'Sheet',sh{1},'WriteRowNames',true);
   writetable(tbl_stats,filename,'Sheet',sh{2});
   if ispc
        try
            objExcel = actxserver('Excel.Application');
            objExcel.Workbooks.Open(filename); 
            % MRZ: delte default sheets in freshly created .xls
            [~, sheets] = xlsfinfo(filename);
            sheetNames2remove = setdiff(sheets,sh); 
            for i = 1:numel(sheetNames2remove)
                objExcel.ActiveWorkbook.Worksheets.Item(sheetNames2remove{i}).Delete;
            end
            objExcel.ActiveWorkbook.Save
            objExcel.ActiveWorkbook.Close(filename);
            delete(objExcel);
        catch
        end         
    end
end
