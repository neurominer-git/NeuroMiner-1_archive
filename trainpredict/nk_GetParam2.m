function [param, model] = nk_GetParam2(Y, label, Params, ModelOnly)
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
[param, model] = feval( TRAINFUNC, Y, label, ModelOnly, Params );

