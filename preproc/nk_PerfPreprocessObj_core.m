% =========================================================================
% FORMAT [InputParam, TrainedParam] = nk_PerfPreprocess_core(InputParam, 
%                                                             TrParam, act)
% =========================================================================
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 1/2016

function [SrcParam, InputParam, TrParam] = nk_PerfPreprocessObj_core(SrcParam, InputParam, TrParam, act)

global VERBOSE TEMPL

%% Prepare and check input params
paramfl = false; tsfl = false; trfl = false; cfl = false;

% If training params are empty they have to be computed right now
if ~exist('TrParam','var') || isempty(TrParam),
    if VERBOSE; fprintf('\nNo Training parameters found. Training mode enabled.'); end
    TrParam = []; 
else
    fprintf('\nTraining parameters found. Out-of-sample mode enabled.');   
    paramfl = true; 
end
% Is data there to be processed?
if ~exist('InputParam','var') || isempty(InputParam), return; end

% Training data?
if isfield(InputParam,'Tr'), trfl=true; end

% Test data and multiple test data sets?
nTs = 0; if isfield(InputParam,'Ts'), tsfl=true; if iscell(InputParam.Ts), nTs = size(InputParam.Ts,2); end; end
nC = 0; if isfield(InputParam,'C'), cfl=true; if iscell(InputParam.C), nC = size(InputParam.C,2); end; end

if ~trfl && ~tsfl, error('\nNo Preprocessing performed because training and test data are missing.\nCheck your parameters!!!'); end

% Is action sequence provided? 
if iscell(act), nact = length(act); else nact = 1; acti = act; end

% Do we training parameters and does training param sequence match action sequence 
if iscell(TrParam), nTrParam = length(TrParam); else, nTrParam = 1; TrParami = []; end

if nact ~= nTrParam && paramfl, if VERBOSE;fprintf('\nLength of action sequence has to match training parameters. Abort!!!'); end; return; end
sumElapsed=0;

if nact>1, fprintf('\nCV1 [%g, %g]: Start preprocessing sequence.', InputParam.CV1perm, InputParam.CV1fold); end

for i=1:nact
    
    ActParam = struct('trfl', trfl, 'paramfl', paramfl, 'tsfl', tsfl, 'cfl', cfl, 'nTs', nTs', 'nC', nC,'i',i);
        
    if iscell(act), acti = act{i}; end
    % Out-of-sample mode: get current training parameters
    if iscell(TrParam) && paramfl, TrParami = TrParam{i}; end
    
    tic; if VERBOSE;fprintf('\nStep #%g/%g:',i,nact); end
    funcstr = sprintf('act_%s', acti); 
    % Check whether multiple shelfs of data are available
    if iscell(InputParam.Tr) 
        nP = numel(InputParam.Tr);
        if isfield(InputParam.P{i},'opt') && ~isempty(InputParam.P{i}.opt)
            O = InputParam.P{i}.opt;
            nO = size(InputParam.P{i}.opt,1);
        else
            O = [];
            nO = 1;
        end
        nComb = nP * nO; 
        Tr = cell(nComb,1); 
        if isfield(InputParam,'Yw'), Yw = cell(nComb,1); end
        if nTs>0, Ts = cell(nComb,nTs); end
        TrParamij = cell(nComb,1); ll = 1;
        % Prepare Input parameter container for multi-parameter operations
        InputParamj.P               = InputParam.P;
        InputParamj.CV1perm         = InputParam.CV1perm;
        InputParamj.CV1fold         = InputParam.CV1fold;
        InputParamj.curclass        = InputParam.curclass;
        if nTs > 1, InputParamj.Ts  = cell(1,nTs);   end
        % Loop through parameter space
        for j=1:nP
            ActParam.j = j;
            % Get training and test data from current parameter shelf
            if numel(InputParam.Tr) == 1
               InputParamj.Tr = InputParam.Tr{1}; 
            else
               InputParamj.Tr = InputParam.Tr{j}; 
            end
            if isfield(InputParam,'Yw'), 
                try
                    InputParamj.Yw = InputParam.Yw{j}; 
                catch
                    InputParamj.Yw = InputParam.Yw;
                end
            end
            if ~isempty(TEMPL), 
                if isstruct(TEMPL.Param{InputParam.curclass}{i})
                    ActParam.Templ = TEMPL.Param{InputParam.curclass}{i};
                else
                    ActParam.Templ = TEMPL.Param{InputParam.curclass}{i}{j}; 
                end
            end
            if nTs>1
                for k = 1:nTs
                    if size( InputParam.Ts,1)==1
                        InputParamj.Ts{k} = InputParam.Ts{1,k};
                    else
                        InputParamj.Ts{k} = InputParam.Ts{j,k};
                    end
                end
            elseif nTs==1
                if size( InputParam.Ts,1)==1
                    InputParamj.Ts = InputParam.Ts{1};
                else
                    InputParamj.Ts = InputParam.Ts{j};
                end
            end

            % ... and now apply multiple parameters to shelf thus creating
            % multi-shelf versions of the data
            for l = 1: nO
                if ~isempty(O)
                    ActParam.opt = O(l,:);
                elseif isfield(ActParam,'opt'), 
                    ActParam = rmfield(ActParam,'opt');
                end
                if paramfl, llTrParam = TrParami{ll}; else, llTrParam = []; end
                [ SrcParam, Out, TrParamij{ll}, ActParam ] = feval( funcstr, SrcParam, InputParamj, TrParam, llTrParam, ActParam );
%                 if isfield(TrParamij{ll},'indNonRem')
%                     ActParam.indNonRem{i,j,l} = TrParamij{ll}.indNonRem;
%                 end
%                 if ~isempty(O),
%                     TrParamij{ll}.i_opt = O(l,:);
%                 end
                Tr{ll} = Out.Tr;
                if isfield(Out,'Yw'), Yw{ll} = Out.Yw; end
                if nTs>1
                    for k = 1:nTs,
                        Ts{ll,k} = Out.Ts{k};
                    end
                elseif nTs==1
                    Ts{ll} = Out.Ts;
                end
                ll=ll+1;
            end
        end
        InputParam.Tr = Tr; 
        if isfield(InputParam,'Yw'), InputParam.Yw = Yw; end
        if nTs > 0; InputParam.Ts = Ts; end
        TrParami = TrParamij;
        
    elseif isfield(InputParam.P{i},'opt') && ~isempty(InputParam.P{i}.opt)
        O = InputParam.P{i}.opt;
        nO = size(InputParam.P{i}.opt,1);
        Tr = cell(nO,1); 
        if isfield(InputParam,'Yw'), Yw = cell(nO,1); end
        if nTs > 0, Ts = cell(nO,nTs); end
        TrParamij = cell(nO,1); ll = 1;
        for l = 1: nO
            ActParam.opt = O(l,:);
            if ~isempty(TEMPL), ActParam.Templ = TEMPL.Param{InputParam.curclass}{i}; end
            if paramfl, llTrParam = TrParami{ll}; else, llTrParam = []; end
            [ SrcParam, Out, TrParamij{ll}, ActParam ] = feval( funcstr, SrcParam, InputParam, TrParam, llTrParam, ActParam );
            TrParamij{ll}.i_opt = O(l,:);
            Tr{ll} = Out.Tr;
            if isfield(Out,'Yw'), Yw{ll} = Out.Yw; end
            if nTs>1
                for k = 1:nTs
                    Ts{ll,k} = Out.Ts{k};
                end
            elseif nTs==1
                Ts{ll} = Out.Ts;
            end
            ll=ll+1;
        end
        InputParam.Tr = Tr; 
        if isfield(InputParam,'Yw'), InputParam.Yw = Yw; end
        if nTs > 0; InputParam.Ts = Ts; end
        TrParami = TrParamij;
    else
        if ~isempty(TEMPL), ActParam.Templ = TEMPL.Param{InputParam.curclass}{i}; end
        [ SrcParam, InputParam, TrParami ] = feval( funcstr, SrcParam, InputParam, TrParam, TrParami, ActParam );
    end
    
    tElapsed=toc; if VERBOSE, fprintf(' Done in %1.2fs.',tElapsed); end; sumElapsed = sumElapsed + tElapsed;
    
    % No Out-of-sample mode => build up TrParam structure array
    if ~paramfl, TrParam{i} = TrParami; end; TrParami = []; 
    
end
if nact > 1, fprintf('\tCompleted in %1.2fs.',sumElapsed); end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam] = act_impute(SrcParam, InputParam, ~, TrParami, actparam)
%global VERBOSE
trfl     = actparam.trfl;
tsfl     = actparam.tsfl;
%paramfl  = actparam.paramfl;
i        = actparam.i;
tsproc   = false;
InputParam.P{i}.IMPUTE.X = InputParam.Tr;
if trfl, 
    [InputParam.Tr, TrParami] = nk_PerfImputeObj(InputParam.Tr, InputParam.P{i}.IMPUTE); 
    if tsfl, tsproc = true; end  
end
if tsproc, InputParam.Ts = nk_PerfImputeObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [ SrcParam, InputParam, TrParami, actparam] = act_labelimpute(SrcParam, InputParam, ~, ~, actparam)
global MODEFL
i  = actparam.i;
InputParam.P{i}.LABELIMPUTE.X = InputParam.Tr;
switch MODEFL
    case 'classification'
        L = SrcParam.MultiTrainLabel;
    case 'regression'
        L = SrcParam.TrainLabel;
end
[SrcParam.TrL_imputed, TrParami] = nk_PerfImputeLabelObj(L, InputParam.P{i}.LABELIMPUTE); 

end

% =========================================================================
function [ SrcParam, InputParam, TrParami, actparam ] = act_scale(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;
if ~isfield(InputParam.P{i}.SCALE,'AcMatFl'), InputParam.P{i}.SCALE.AcMatFl = false; end
if ~isfield(InputParam.P{i}.SCALE,'ZeroOne'), InputParam.P{i}.SCALE.ZeroOne = 1;     end
if ~isfield(InputParam.P{i}.SCALE,'zerooutflag'), InputParam.P{i}.SCALE.zerooutflag = 1;     end
if VERBOSE
    if InputParam.P{i}.SCALE.ZeroOne == 1, rangestr = '[0, 1]'; else rangestr = '[-1, 1]'; end
    if InputParam.P{i}.SCALE.AcMatFl, fprintf('\tScaling whole matrix to %s ...', rangestr); else fprintf('\tScaling each feature to %s ...', rangestr); end
end
if paramfl && tsfl && isfield(TrParami,'minY') && ...
        isfield(TrParami,'maxY'), tsproc = true;
else
    if trfl, [InputParam.Tr, TrParami] = nk_PerfScaleObj(InputParam.Tr, InputParam.P{i}.SCALE); end
    if tsfl, tsproc = true; end  
end
if tsproc, InputParam.Ts = nk_PerfScaleObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [ SrcParam, InputParam, TrParami, actparam ] = act_normalize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc = false;

if VERBOSE; fprintf('\tNormalizing to group mean(s) ...'); end
            
if paramfl && tsfl && isfield(TrParami,'meanY') && isfield(InputParam.P{i},'TsInd')
    tsproc = true;
elseif trfl, 
    [InputParam.Tr, TrParami] = nk_PerfNormObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl, tsproc = true; end 
else
    if VERBOSE; fprintf('\n not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfNormObj(InputParam.Ts, TrParami); end

end
% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_copyfeat(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

if VERBOSE;fprintf('\tCopy features ...'); end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_unitnormalize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if VERBOSE; fprintf('\tNormalizing to unit vector ...'); end 

if paramfl && tsfl && isfield(TrParami,'normY') 
    tsproc = true;
elseif isfield(InputParam.P{i},'METHOD')
    if trfl, [InputParam.Tr, TrParami] = nk_PerfUnitNormObj(InputParam.Tr, InputParam.P{i}); end
    if tsfl, tsproc = true; end 
else
    if VERBOSE;fprintf('\n not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfUnitNormObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_remmeandiff(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if VERBOSE; fprintf('\tCorrect groups for offsets from global mean(s) ...'); end

if paramfl && tsfl && isfield(TrParami,'meanY') && isfield(TrParami,'meanG')
    tsproc = true;
elseif trfl, 
    [InputParam.Tr, TrParami] = nk_PerfRemMeanDiffObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl, tsproc = true; end 
else
    if VERBOSE;fprintf(' not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfRemMeanDiffObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_discretize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if VERBOSE;
    fprintf('\tDiscretizing to range [ %g : %g : %g ] ...', ...
        InputParam.P{i}.DISCRET.binstart, ...
        InputParam.P{i}.DISCRET.binsteps, ...
        InputParam.P{i}.DISCRET.binstop)
end
if paramfl && tsfl && isfield(TrParami,'mY') && isfield(TrParami,'sY'), 
    tsproc = true;
elseif trfl,
    [InputParam.Tr, TrParami] = nk_PerfDiscretizeObj(InputParam.Tr, InputParam.P{i});
    if tsfl, tsproc = true; end
else
    if VERBOSE;fprintf(' not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfDiscretizeObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_symbolize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
nTs     = actparam.nTs;
tsproc = false;    

if VERBOSE; fprintf('\tSymbolizing ...'); end

if paramfl && tsfl && isfield(TrParami,'sBin')
    tsproc=true;
elseif trfl, 
    [ InputParam.Tr, TrParami ] = nk_PerfDiscretizeObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl, tsproc = true; end;
else
    if VERBOSE;fprintf(' not performed.'); end
end
if tsproc, InputParam.Ts = nk_PerfDiscretizeObj(InputParam.Ts, TrParami); end

end   

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_standardize(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc = false;            

if isfield(actparam,'opt')
    InputParam.P{i}.WINSOPT = actparam.opt;
end

if VERBOSE; fprintf('\tStandardizing data ...'); end
if paramfl && tsfl && isfield(TrParami,'meanY') && isfield(TrParami,'stdY')
    tsproc = true;
else
    if trfl, [InputParam.Tr, TrParami] = nk_PerfStandardizeObj(InputParam.Tr, InputParam.P{i}); end
    if tsfl, tsproc = true; end
end
if tsproc, InputParam.Ts = nk_PerfStandardizeObj(InputParam.Ts, TrParami); end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_reducedim(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'DIMOPT') && ~isempty(actparam.DIMOPT)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.DR.opt = actparam.DIMOPT;
elseif isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.DR.opt = actparam.opt;
end 

% We need the template parameters for the optional Procrustes rotation
if isfield(actparam,'Templ') 
    if isfield(actparam.Templ,'mpp')
        InputParam.P{i}.DR.mpp_template = actparam.Templ.mpp;
    elseif iscell(actparam.Templ) && isfield(actparam.Templ{1},'mpp')
        InputParam.P{i}.DR.mpp_template = actparam.Templ{1}.mpp;
    end
end

InputParam.P{i}.TrX = InputParam.Tr;

if VERBOSE    
    if trfl, 
        featnum = size(InputParam.Tr,2); 
    elseif tsfl
        if iscell(InputParam.Ts)
            featnum = size(InputParam.Ts{1},2); 
        else
            featnum = size(InputParam.Ts,2); 
        end
    end
    fprintf('\tDimensionality reduction [%s]: %g input features -> ', ...
                                    InputParam.P{i}.DR.RedMode, featnum); 
end

if paramfl && tsfl && isfield(InputParam.P{i},'DR') && ...
        isfield(TrParami,'mpp')
    tsproc = true;
elseif trfl, 
    [InputParam.Tr, TrParami] = nk_PerfRedObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl,tsproc = true; end   
end 

if tsproc, InputParam.Ts = nk_PerfRedObj(InputParam.Ts, TrParami); end

end

function [SrcParam,InputParam, TrParami, actparam ] = act_extdim(SrcParam, InputParam, TrParam, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.opt = actparam.opt;
end 

% Extract weight vector from last TrParam
if isfield(actparam,'j'); j = actparam.j; end
if exist('j','var')
    if isfield(TrParam{i-1}{j},'mpp')
        InputParam.P{i}.mpp = TrParam{i-1}{j}.mpp;
    end
else
    if isfield(TrParam{i-1},'mpp')
        InputParam.P{i}.mpp = TrParam{i-1}.mpp;
    end
end

if VERBOSE    
    if trfl, 
        featnum = size(InputParam.Tr,2); 
    elseif tsfl
        if iscell(InputParam.Ts)
            featnum = size(InputParam.Ts{1},2); 
        else
            featnum = size(InputParam.Ts,2); 
        end
    end
    fprintf('\tExtraction of subspaces from reduced data: %g input features -> ', featnum); 
end

if paramfl && tsfl && isfield(TrParami,'indNonRem')
    tsproc = true;
elseif trfl, 
    [InputParam.Tr, TrParami] = nk_PerfExtDimObj(InputParam.Tr, InputParam.P{i}); 
    if tsfl,tsproc = true; end   
end 

if tsproc, InputParam.Ts = nk_PerfExtDimObj(InputParam.Ts, TrParami); end

end
% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_correctnuis(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;            

IN = InputParam.P{i}; 

IN.revertflag = 0; if isfield(InputParam.P{i},'COVDIR') && InputParam.P{i}.COVDIR ==1, IN.revertflag = 1; end
IN.nointercept = 1; if isfield(InputParam.P{i},'INTERCEPT'), IN.nointercept = ~InputParam.P{i}.INTERCEPT; end

if VERBOSE
    fprintf('\tPartial correlations ...'); 
    if IN.revertflag, fprintf(' introduce covariate effects ...'); else fprintf(' remove covariate effects ...'); end
    if IN.nointercept, fprintf(' intecept disabled ...'); else fprintf(' intercept enabled ...'); end
end

if isfield(InputParam.P{i},'BETAEXT')
    if VERBOSE;fprintf(' external beta found ...'); end 
    IN.beta = InputParam.P{i}.BETAEXT; 
else
    if isfield(InputParam.P{i},'SUBGROUP')
        IN.subgroup = InputParam.P{i}.SUBGROUP;
        if VERBOSE;fprintf(' subgroup beta computation ...'); end
    end
end
if paramfl && tsfl && isfield(TrParami,'beta') 
    % Out-of-sample mode, used TsCovars and stored beta
    tsproc = true;
elseif trfl, 
    [InputParam.Tr, TrParami] = nk_PartialCorrelationsObj(InputParam.Tr, IN);
    if tsfl, tsproc = true; end
else
    if VERBOSE;fprintf('not performed'); end
end

if tsproc, InputParam.Ts = nk_PartialCorrelationsObj(InputParam.Ts, TrParami); end 
end

% =========================================================================
function [SrcParam, InputParam, TrParami, actparam ] = act_elimzero(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;  

if isfield(actparam,'opt')
    InputParam.P{i}.PRUNE.perc = actparam.opt;
end

if paramfl && tsfl && isfield(TrParami,'NonPruneVec') 
    tsproc = true;
elseif trfl, 
    if VERBOSE;fprintf('\tAttribute pruning ...'); end
    [InputParam.Tr, TrParami] = nk_PerfElimZeroObj(InputParam.Tr, InputParam.P{i}.PRUNE);
    % All subsequent processing steps that use fixed column indices have to
    % be adjusted to the pruned matrix 
    for z=i+1:numel(InputParam.P)
        if isfield(InputParam.P{z},'IMPUTE') && ~isempty(InputParam.P{z}.IMPUTE.blockind)
            InputParam.P{z}.IMPUTE.blockind = InputParam.P{z}.IMPUTE.blockind(TrParami.NonPruneVec);
        end
    end
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = nk_PerfElimZeroObj(InputParam.Ts, TrParami); end
end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_rankfeat(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE

i  = actparam.i;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.RANK.opt = actparam.opt;
end
if isfield(InputParam,'Yw') && strcmp(InputParam.P{i}.RANK.algostr,'extern')
   InputParam.P{i}.RANK.EXTERN = InputParam.Yw;
end

if ~isfield(TrParami,'W') || isempty(TrParami.W)
    if VERBOSE;fprintf('\tCompute Feature Weighting ...'); end
    TrParami = nk_PerfFeatRankObj(InputParam.Tr, InputParam.P{i}.RANK);
    actparam.RANK.W = TrParami.W;
end

end
% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_extfeat(SrcParam, InputParam, TrParam, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.W_ACT.opt = actparam.opt;
end

% Extract weight vector from last TrParam
if isfield(actparam,'j'); j = actparam.j; end
InputParam.P{i}.plsdev = false; 

if exist('j','var')
    InputParam.P{i}.W = TrParam{i-1}{j}.W;
else
    InputParam.P{i}.W = TrParam{i-1}.W;
end

if VERBOSE,fprintf('\tRanking-based feature extraction ...'); end
if paramfl && tsfl 
    % Out-of-sample mode, used stored weights 
    tsproc = true;
elseif trfl, 
   % Training mode, learn beta and optionally apply it to
    % out-of-sample data
    [ InputParam.Tr, TrParami ] = nk_PerfWActObj(InputParam.Tr, InputParam.P{i});
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = nk_PerfWActObj(InputParam.Ts, TrParami ); end

if isfield(InputParam,'Yw') && isfield(TrParami,'ind') && isnumeric(InputParam.Yw)
    InputParam.Yw = InputParam.Yw(TrParami.ind);
end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_remvarcomp(SrcParam, InputParam, ~, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.REMVARCOMP.corrthresh = nk_ReturnParam('CorrThresh', InputParam.P{i}.REMVARCOMP.Params_desc, actparam.opt); 
    InputParam.P{i}.REMVARCOMP.DR.dims = nk_ReturnParam('Dims', InputParam.P{i}.REMVARCOMP.Params_desc, actparam.opt); 
end

% We need the template parameters for the optional Procrustes rotation
if isfield(actparam,'Templ') && isfield(actparam.Templ,'mpp')
   InputParam.P{i}.REMVARCOMP.DR.mpp_template = actparam.Templ.mpp;
end

if VERBOSE,fprintf('\tExtract variance components ...'); end
if paramfl && tsfl 
    % Out-of-sample mode, used stored weights 
    tsproc = true;
elseif trfl, 
   % Training mode, learn beta and optionally apply it to
    % out-of-sample data
    [ InputParam.Tr, TrParami ] = nk_PerfAdjForCovarsUsingPCAObj(InputParam.Tr, InputParam.P{i}.REMVARCOMP, InputParam.Tr);
    
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = nk_PerfAdjForCovarsUsingPCAObj(InputParam.Ts, TrParami ); end

end

% =========================================================================
function [SrcParam,InputParam, TrParami, actparam ] = act_devmap(SrcParam, InputParam, TrParam, TrParami, actparam)
global VERBOSE
trfl    = actparam.trfl;
tsfl    = actparam.tsfl;
paramfl = actparam.paramfl;
i       = actparam.i;
tsproc  = false;

if isfield(actparam,'opt') && ~isempty(actparam.opt)
    % Overrun user-defined dimensionality!
    InputParam.P{i}.DEVMAP.cu = nk_ReturnParam('spls-cu', InputParam.P{i}.DEVMAP.Params_desc, actparam.opt); 
    InputParam.P{i}.DEVMAP.cv = nk_ReturnParam('spls-cv', InputParam.P{i}.DEVMAP.Params_desc, actparam.opt); 
end

%InputParam.P{i}.DEVMAP.loo = true;

if VERBOSE,fprintf('\tDeviation-based feature weighting ...'); end
if paramfl && tsfl 
    % Out-of-sample mode, used stored weights 
    tsproc = true;
elseif trfl, 
   % Training mode, learn beta and optionally apply it to
    % out-of-sample data
    [ InputParam.Tr, TrParami ] = nk_PerfDevMapObj(InputParam.Tr, InputParam.P{i});
    if tsfl, tsproc = true; end
end

if tsproc, InputParam.Ts = nk_PerfDevMapObj(InputParam.Ts, TrParami ); end

end
