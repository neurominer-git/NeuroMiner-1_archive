function [sY, IN] = nk_PartialCorrelationsObj( Y, IN )
% =========================================================================
% FORMAT [Y, IN] = nk_PartialCorrelationsObj(Y, IN)
% =========================================================================
% Remove nuisance effects IN.G from Y (optionally, using a predefined beta)
%
% I\O Arguments:
% -------------------------------------------------------------------------
% Y                 : M cases x N features data matrix
% IN.G              : The covariate(s) to regress out from Y
% IN.nointercept    : Include an intercept in the model or not
% IN.subgroup       : Index vector of cases in Y to compute the beta(s) from
% IN.beta           : The estimated beta coefficients
% IN.revertflag     : Increase (=1) or attenuate (=0) IN.G effects 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10 / 2015

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y));
    for i=1:numel(Y), 
        if isfield(IN,'beta') && ~isempty(IN.beta)
            IN.G = IN.TsCovars{i};
        else
            IN.G = IN.TrCovars{i};
        end
        [sY{i}, IN] = PartialCorrelationsObj(Y{i}, IN); 
    end
else
    if isfield(IN,'beta') && ~isempty(IN.beta)
        if iscell(IN.TsCovars)
            IN.G = IN.TrCovars;
        else
            IN.G = IN.TsCovars;
        end
    else
        IN.G = IN.TrCovars;
    end
    [ sY, IN ] = PartialCorrelationsObj( Y, IN );
end
% =========================================================================

function [Y, IN] = PartialCorrelationsObj( Y, IN )

if isempty(IN),eIN=true; else eIN=false; end

if eIN|| ~isfield(IN,'G') || isempty(IN.G), error('No covariates defined in parameter structure'), end

if eIN || (~isfield(IN,'nointercept') || isempty(IN.nointercept) || ~IN.nointercept ) 
     %Create intercept vecotr
    interceptflag = true;
    intercept = ones(size(IN.G,1),1);
    %Check if intercept is already included in matrix to avoid double
    %intercept removal
    if isfield(IN,'beta') && ~isempty(IN.beta)
        if size(IN.beta,1) == size(IN.G,2), interceptflag = false; end
    end
else
    interceptflag = false;
end

if interceptflag
    %fprintf(' ... adding intercept to covariate matrix')
    IN.G = [intercept IN.G];
end

if eIN || ~isfield(IN,'beta') || isempty(IN.beta), 
    if ~isfield(IN,'subgroup') || isempty(IN.subgroup)
        % Compute IN.beta from entire population
        IN.beta = pinv(IN.G) * Y; 
    else
        % Compute IN.beta from a subgroup of observations
        IN.beta = pinv(IN.G(IN.subgroup,:)) * Y(IN.subgroup,:);
    end
end
if eIN || ~isfield(IN,'revertflag') || isempty(IN.revertflag) || ~IN.revertflag
    try 
        Y = Y - IN.G * IN.beta;
    catch
        fprintf('problem')
    end
else
    Y = Y + IN.G * IN.beta;
end

return
