% =========================================================================
% FORMAT param = CC(expected, predicted)
% =========================================================================
% Compute Correlation Coefficient
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 09/2011

function param = CC(expected, predicted)
if isempty(expected), param = []; return; end
param = corrcoef(expected,predicted);
    param = param(2);
end