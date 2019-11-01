% =========================================================================
% FORMAT function [param, model] = nk_GetParam_LIBSVM(Y, label, SlackParam, ...
%                                                KernParam, ModelOnly)
% =========================================================================
% Train LIBSVM models and evaluate their performance using Y & label, 
% SlackParam, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2012

function [param, model] = nk_GetParam_LIBSVM(Y, label, SlackParam, ...
                                                KernParam, ModelOnly)
                                            
global SVM CMDSTR EVALFUNC LIBSVMTRAIN MODEFL                         

param = [];flw = 0;
if strcmp(LIBSVMTRAIN,'svmtrain312'), flw = 1; end

cmdstr = [CMDSTR.simplemodel SlackParam CMDSTR.KernStr KernParam ' -q'];

% Check if weighting is necessary
if SVM.LIBSVM.Weighting 
    flg = 0;
    switch MODEFL
        case 'classification'
            % Works only for C-SVC (not nu-SVC) as +C / -C slacks are treated separately
            npos = sum(label == 1); nneg = numel(label) - npos; bal = npos / nneg;
            if  abs(bal-1) > 0.1
                fac2 = 1; fac1 = 1/bal;
                cmdstr = [ cmdstr ' -w1 ' num2str(fac1) ' -w-1 ' num2str(fac2)]; 
            end
        case 'regression'
            W = nk_WeigthDataInstanceHisto(label);
    end
elseif flw
    W = ones(numel(label),1);
end

if iscell(Y) 
   
    % MKL-based learning not implemented yet
   
else % Univariate case
    if size(label,1) ~= size(Y,1), label = label'; end
    if flw
         model = feval( LIBSVMTRAIN, W, label, Y, cmdstr );
    else
        model = feval( LIBSVMTRAIN, label, Y, cmdstr );
    end
    if ~ModelOnly
        [param.target, param.dec_values] = nk_GetTestPerf_LIBSVM([], Y, label, model) ;
        param.val = feval(EVALFUNC, label, param.target);
    end
end