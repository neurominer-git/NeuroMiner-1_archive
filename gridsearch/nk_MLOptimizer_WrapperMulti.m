function [R, optmodel] = nk_MLOptimizer_WrapperMulti(Y, mY, L, bL, MultiL, mYnew, Lnew, MultiLnew, ngroups, Ps, FullParam, SubFeat)

global RFE SVM

if nargin < 12, SubFeat = true(1,size(Y,2)); end
ActStr = {'Tr', 'CV', 'TrCV'};

% Remove cases which are completely NaN
for i=1:numel(Y)
    [Y{i}, L{i}] = nk_ManageNanCases(Y{i}, L{i});
    
    % Run ADASYN if needed
    if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
        [Y{i}, L{i}] = nk_PerfADASYN(Y{i}, L{i}, SVM.ADASYN);
    end
    
    if i==1
        [mY{i}, MultiL, I] = nk_ManageNanCases(mY{i}, MultiL); 
        [mYnew{i}, Lnew{i}, Inew] = nk_ManageNanCases(mYnew{i}, Lnew{i}); 
        MultiLnew(Inew)=[]; bL{i}(I)=[];
    else
        mY{i}(I,:)=[]; bL{i}(I)=[];
        mYnew{i}(Inew,:)=[]; Lnew{i}(Inew)=[]; 
    end
    
end

switch RFE.Wrapper.type
    %% GREEDY FORWARD/BACKWARD FEATURE SEARCH
    case 1 
        switch RFE.Wrapper.GreedySearch.Direction
            case 1
                % Forward selection using argmax => CV1 => test data performance as criterion 
                rfefun = 'rfe_forward_multi';                
            %%% BACKWARD ELIMINATION %%%
            case 2 % Backward elimination using argmax => CV1 => test &train data performance as criterion 
                rfefun = 'rfe_backward_multi';
        end
        [optparam, optind, optfound, optmodel] = feval(rfefun, Y, mY, L, bL, MultiL, mYnew, Lnew, MultiLnew, Ps, SubFeat, FullParam, ngroups, ActStr{RFE.Wrapper.datamode});
    
    %%% SIMULATED ANNEALING %%%
    case 2
        [optparam, optind, optfound, optmodel] = ...
            nk_SimAnneal(Y, label, Ynew, labelnew, Ps, SubFeat, FullParam, ActStr{RFE.Wrapper.datamode});
end
% Transfer params to output structure
R.found                   = optfound;
R.FeatureIndex            = optind;
if nargin == 8,
    R.SubFeatIndex        = false(size(SubFeat));
    R.SubFeatIndex(optind)= SubFeat(optind); 
end
R.OptimParam              = optparam;


