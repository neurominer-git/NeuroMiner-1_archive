% ==========================================================================
% FORMAT [param, model] = nk_GetParam_LIBLIN(Y, label, SlackParam, ~, ...
%                                           ModelOnly)
% ==========================================================================
% Train LIBLINEAR models and evaluate their performance using Y & label, 
% SlackParam,
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2012

function [param, model] = nk_GetParam_LIBLIN(Y, label, SlackParam, ~, ModelOnly)
global SVM CMDSTR MODEFL EVALFUNC                          

param = [];
cmd = [ '-c ' SlackParam ];
       
% Check if weighting is necessary
cmd = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmd); cmd = [ cmd CMDSTR.simplemodel ];

if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    sY=sparse(Y);
    model = train_liblin(label, sY, cmd);
    if ~ModelOnly
        [param.target, ...
            param.val, ...
            param.dec_values] = predict_liblin(label, sY, model, SVM.LIBLIN.b);
        param.val = feval(EVALFUNC, label, param.target);
    end
end