function [SVMfeat, optmodel] = nk_MLOptimizer_Wrapper(Y, label, Ynew, labelnew, Ps, FullParam, SubFeat)

global RFE
if nargin < 8, SubFeat = true(1,size(Y,2)); end
ActStr = {'Tr', 'CV', 'TrCV'};
% Remove cases which are completely NaN
[Y, label] = nk_ManageNanCases(Y, label);
[Ynew, labelnew] = nk_ManageNanCases(Ynew, labelnew);
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
SVMfeat.found                   = optfound;
SVMfeat.FeatureIndex            = optind;
if nargin == 8,
    SVMfeat.SubFeatIndex        = false(size(SubFeat));
    SVMfeat.SubFeatIndex(optind)= SubFeat(optind); 
end
SVMfeat.OptimParam              = optparam;


