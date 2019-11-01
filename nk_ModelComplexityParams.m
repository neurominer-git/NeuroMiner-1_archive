function [sumNSV, meanNSV, stdNSV, rSV] = nk_ModelComplexityParams(M, S)
% function [sumNSV, meanNSV, stdNSV] = nk_ModelComplexityParams(M)
%
% Compute model complexity parameters, including ...:
% sum of support vectors (SV):      sumNSV
% mean number of SV:                meanNSV
% std number of SV:                 stdNSV
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 12/2009

global SVM

[ix,jx,nclass] = size(M); 

sumNSV = zeros(nclass,1);
meanNSV = zeros(nclass,1);
stdNSV = zeros(nclass,1);
rSV = zeros(nclass,1);

for curclass=1:nclass
    
    NSV=[];
    switch SVM.prog

        case {'LIBSVM','MikRVM','MVTRVR','MSTOOL'}
            for i=1:ix
                for j=1:jx
                    for l=1:length(M{i,j,curclass})
                       
                        NSV = [NSV; M{i,j,curclass}{l}.totalSV];
                        
                    end
                end
            end
        case 'MKLRVM'
            for i=1:ix
                for j=1:jx
                    for l=1:length(M{i,j,curclass})
                       
                        NSV = [NSV; M{i,j,curclass}{l}.N_prototypical];
                        
                    end
                end
            end
        case 'SEQOPT'
            for i=1:ix
                for j=1:jx
                    for l=1:length(M{i,j,curclass})                       
                        NSV = [NSV; mean(M{i,j,curclass}{l}.examsfreq(end))];
                    end
                end
            end
        otherwise
            NSV = S(curclass);
    end

    sumNSV(curclass)  = sum(NSV);
    meanNSV(curclass) = mean(NSV);
    stdNSV(curclass)  = std(NSV);
    switch SVM.prog
        case 'SEQOPT'
            rSV(curclass)     = meanNSV(curclass);
        otherwise
            rSV(curclass)     = meanNSV(curclass)/S(curclass)*100;
    end
end

