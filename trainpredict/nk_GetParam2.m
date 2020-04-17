function [param, model] = nk_GetParam2(Y, label, Params, ModelOnly)
% =========================================================================
% FORMAT [param, model] = nk_GetParam2(Y, label, Params, ModelOnly)
% =========================================================================
% Generic interface function to the algorithm-specific training modules
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2017
global TRAINFUNC

% Remove cases which are completely NaN
[Y, label] = nk_ManageNanCases(Y, label);

% Pass training matrix and labels to training module
[param, model] = feval( TRAINFUNC, Y, label, ModelOnly, Params );

