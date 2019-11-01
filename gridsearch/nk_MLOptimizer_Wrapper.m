function [R, optmodel] = nk_MLOptimizer_Wrapper(Y, label, Ynew, labelnew, Ps, FullParam, SubFeat)

global RFE SVM

if nargin < 8, SubFeat = true(1,size(Y,2)); end
ActStr = {'Tr', 'CV', 'TrCV'};

% Remove cases which are completely NaN
[Y, label] = nk_ManageNanCases(Y, label);
[Ynew, labelnew] = nk_ManageNanCases(Ynew, labelnew);

% Run ADASYN if needed
if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
    [Y, label] = nk_PerfADASYN(Y, label , SVM.ADASYN);
end

switch RFE.Wrapper.type
    %% GREEDY FORWARD/BACKWARD FEATURE SEARCH
    case 1 
        switch RFE.Wrapper.GreedySearch.Direction
            case 1
                % Forward selection using argmax => CV1 => test data performance as criterion 
                rfefun = 'rfe_forward';                
            %%% BACKWARD ELIMINATION %%%
            case 2 % Backward elimination using argmax => CV1 => test &train data performance as criterion 
                rfefun = 'rfe_backward';
        end
        [optparam, optind, optfound, optmodel] = feval(rfefun, Y, label, Ynew, labelnew, Ps, SubFeat, FullParam, ActStr{RFE.Wrapper.datamode});
    
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


