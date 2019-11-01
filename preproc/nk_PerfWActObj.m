function [ sY, IN ] = nk_PerfWActObj( Y, IN )
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfWActObj(Y, IN)
% =========================================================================
% 
% Inputs: 
% -------------------------------------------------------------------------
% Y                   : M cases x N features data matrix
% IN                  : Input parameter structure 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2018

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] = PerfWActObj(Y{i}, IN); end
else
    [ sY, IN ] = PerfWActObj(Y, IN );
end

% =========================================================================

function [ Y, IN ] = PerfWActObj( Y, IN )
global VERBOSE

% Defaults
if isempty(IN),          error('Suitable input structure is missing. See the functions'' help for more information.'); end
if ~isfield(IN,'W_ACT'), error('Weighting parameter structure is missing! Add a W_ACT substructure to your input structure.'); end
if ~isfield(IN,'W'),     error('Weighting vector is missing! Add a weighting vector W to your input structure.'); end
if ~isfield(IN,'Mask'), IN.Mask = []; end
if ~isfield(IN,'Thresh') || isempty(IN.Thresh)
    if isfield(IN.W_ACT,'opt')
        Params_desc = IN.W_ACT.Params_desc;
        opt = IN.W_ACT.opt;
    else
        Params_desc=[]; opt =[]; 
    end

     if numel(IN.W_ACT.threshvec)>1 || any(IN.W_ACT.threshvec)
        t = nk_ReturnParam('Thresholds',Params_desc, opt); 
        if ~isempty(t) && sum(any(t)), IN.Thresh = percentile(IN.W, t); end
     else
         IN.Thresh = [];
     end
end

if ~isempty(IN.Thresh),
    if IN.W_ACT.clustflag == 1
        if isempty(IN.Mask)
            Wthresh = zeros(size(IN.W));
            IN.ind = IN.W > IN.Thresh;
            Wthresh(IN.ind) = IN.W(IN.ind);
            IN.Mask = nk_Cluster(Wthresh, IN.W_ACT);
        end
        Y = nk_ExtractClusterData(Y, IN.Mask);
        if VERBOSE, fprintf('\tClusterizing F into %g mean cluster values', nk_Range(IN.Mask)); end
    else
        % Hard feature selection
        IN.ind = IN.W >= IN.Thresh;
        Y = Y(:, IN.ind); 
        if VERBOSE, fprintf('\tSelecting %g / %g feats at %g', size(Y,2), numel(IN.W), IN.Thresh); end
    end
else 
    % Check whether you have to apply an exponential multiplier to W
    ExpMult = nk_ReturnParam('ExpMult',Params_desc, opt); 
    if ExpMult, W = IN.W.^ExpMult; else W = IN.W; end
    % Soft feature selection
    Y = bsxfun(@times, Y, W);
    IN.ind = true(size(W));
    if VERBOSE, fprintf('\tWeighting F'); end
end

