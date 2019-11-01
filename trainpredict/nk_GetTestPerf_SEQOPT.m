% =========================================================================
% FORMAT function [param, model] = nk_GetParam_SEQOPT(Y, label, ModelOnly, 
%                                                                 ...cmdstr)
% =========================================================================
% test sequence predictor model in new data 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2018

function [rs, ds] = nk_GetTestPerf_SEQOPT(~, tXtest, Ytest, model, ~, ~)

global MODEFL

[~, ds ] = nk_OptPredSeq(tXtest, Ytest, model);

switch MODEFL
    case 'classification'
        rs = sign(ds);
    case 'regression'
        rs = ds;
end
