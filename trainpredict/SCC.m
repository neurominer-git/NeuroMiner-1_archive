% =========================================================================
% FORMAT param = SCC(expected, predicted)
% =========================================================================
% Compute Squared Correlation Coefficient, also called Coefficient of
% Determination or alternatively % Explained Variance
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 04/2012

function param = SCC(expected, predicted)
if isempty(expected), param = []; return; end
%param = (1 - ( sum( (expected-predicted).^2 ) / sum( (expected-mean(expected)).^2 ) )) * 100;
param = CC(expected,predicted); param = param ^ 2 *100;
if isnan(param), param=0; end;

end