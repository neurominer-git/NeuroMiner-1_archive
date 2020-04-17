function [P, Perf, d, crit] = nk_MultiDecideDAG(H, Y, Classes, maxfunc, weightflag)
% function [P, Perf, d, crit] = nk_MultiDecideOneVsOne(H, Y, Classes,
% maxfunc, weightflag)
%
% Perform One-Vs-One multi-group classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2011

if iscell(H)
    [iy,jy,nclass]  = size(H);
    d               = cell(iy,jy);
    P               = cell(iy,jy);
    crit            = cell(iy,jy);
    Perf            = zeros(iy,jy);
else
    iy = 1; jy=1;
    nclass          = max(Classes);
end

for k=1:iy % Loop through partitions (in case iscell(H) = true)
    
    for l=1:jy % Loop through folds (in case iscell(H) = true)
        
        if iscell(H)
            nsubj       = size(H{k,l},1);
            M           = OneVsOne(Classes{k,l}, nclass);
            ngroups     = size(M,1);
        else
            nsubj       = size(H,1);
            M           = OneVsOne(Classes, nclass);
            ngroups     = size(M,1);
            d           = zeros(nsubj,ngroups);
        end
        
        Weights = zeros(1,ngroups);
        if iscell(H), tH= []; for curclass=1:nclass, tH = [tH H{k,l,curclass}]; end; else tH = H; end
        
        td = zeros(nsubj,ngroups-1);
        
        for j=1:ngroups-1

                % Get binary classifier for current node (first and last from
                % list)
                mJ          = M(j1,:) == 1 && M(j2,:) == -1;
                ind0        = mJ ~= 0;
                codeword    = repmat(mJ(ind0),nsubj,1);

                % Get number of dichotomizers for the present group
                %Weights(j)  = sum(ind0);

                % Get dichotomizers (e.g binary SVM /RVM) predictions for
                % present codeword
                tHX         = tH(:,ind0);

                % Get index matrix of correct predictions
                indC        = codeword == sign(tHX);

                % Produce some decision value according to the following
                % functions:
                switch maxfunc
                    case 1
                        indE = codeword ~= sign(tHX); % Get index matrix of wrong predictions
                        gC = get_absvals(tHX, indC); gE = get_absvals(tHX, indE);
                    case 2
                        indE = codeword ~= sign(tHX);% Get index matrix of wrong predictions
                        gC = get_absvals(tHX, indC); gE = get_absvals(tHX, indE);
                    case 3
                        indE = codeword ~= sign(tHX);% Get index matrix of wrong predictions
                        gC = get_absvals(tHX, indC); gE = get_absvals(tHX, indE);
                    case 4 
                        gC = sum(indC,2); gE = sum(indE,2);
                end
               if gC > gE


            end

            if weightflag
                [~,mXI] = max(Weights); Weights = Weights./Weights(mXI); td = td./repmat(Weights,nsubj,1);
            end

            if iscell(H)
                d{k,l} = td; [crit{k,l},P{k,l}] = max(d{k,l},[],2);
                errs = P{k,l} ~= Y;
            else
                d = td; [crit,P] = max(d,[],2);
                errs = P ~= Y;
            end

            Perf(k,l) = (1-sum(errs)/length(errs))*100;
        end
    end
end

end

function g = get_absvals(tHX, ind)

g = zeros(size(tHX));
g(ind) = abs(tHX(ind)); % get decision values of correct predictions

end

function oECOC = OneVsOne(Classes, number_classes)

nc = size(Classes,2);
oECOC=zeros([number_classes nc]);
tECOC=zeros([number_classes number_classes*(number_classes-1)/2]);
counter=1;

% Create multi-group coding matrix for one-vs-one
for i=1:number_classes-1
    for j=i+1:number_classes
        tECOC(i,counter)=1;
        tECOC(j,counter)=-1;
        counter=counter+1;
    end
end
% Assign dichotomizers to class vector
for i=1:number_classes
    ind = Classes == i; l = sum(ind);
    oECOC(:,ind) = repmat(tECOC(:,i),1,l);
end

end