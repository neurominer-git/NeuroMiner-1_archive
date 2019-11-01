% =========================================================================
% FORMAT function [param, model] = nk_GetParam_SEQOPT(Y, label, ModelOnly, 
%                                                                 ...cmdstr)
% =========================================================================
% Train LIBSVM models and evaluate their performance using Y & label, 
% SlackParam, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2018

function [param, model] = nk_GetParam2_SEQOPT(Y, label, ModelOnly, PQ)

global EVALFUNC MODEFL SVM

param =[];
model = nk_OptPredSeq(Y, label, [], SVM.SEQOPT.C(PQ.val(1),:), PQ.val(2), [PQ.val(3) PQ.val(4)], EVALFUNC);

if ~ModelOnly
    param.dec_values = model.optD;
    switch MODEFL
        case 'classification'
            param.target = sign(model.optD);
        case 'regression'
            param.target = model.optD;
    end
    param.val = feval(EVALFUNC, label, param.target);
end
