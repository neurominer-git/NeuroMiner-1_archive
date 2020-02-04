function [param, model] = nk_GetParam2(Y, label, Params, ModelOnly, FeatGroups)
% =========================================================================
% FORMAT [param, model] = nk_GetParam2(Y, label, Params, ModelOnly)
% =========================================================================
% Generic interface function to the algorithm-specific training modules
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 11/2018
global TRAINFUNC SVM

%SVM.ADASYN.flag = true;

% Remove cases which are completely NaN
[Y, label] = nk_ManageNanCases(Y, label);

% Run ADASYN if needed
if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
    [Y, label] = nk_PerfADASYN( Y, label , SVM.ADASYN);
end

% Pass training matrix and labels to training module
switch SVM.prog
    case 'SEQOPT'
        if  ~exist('FeatGroups','var') || isempty(FeatGroups)
            [param, model] = feval( TRAINFUNC, Y, label, [], ModelOnly, Params );
        else
            [param, model] = feval( TRAINFUNC, Y, label, FeatGroups, ModelOnly, Params );
        end
    otherwise
        [param, model] = feval( TRAINFUNC, Y, label, ModelOnly, Params );
end

