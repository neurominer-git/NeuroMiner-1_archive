function [act, dat, inp] = nk_VisModelsPrep(act, dat, inp, parentstr)
% =========================================================================
% function [act, dat, inp] = nk_VisModelsPrep(dat, inp, parentstr)
% =========================================================================
% Wrapper function of the NM visualization module, which allows the user to
% interactively chose run-time options for the analysis of model patterns.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 12/2018

global MULTI CV NM

%% Setup defaults
complvec = []; for z=1:numel(dat.analysis), if dat.analysis{z}.status, complvec = [ complvec z ]; end; end

if ~exist('inp','var') || isempty(inp)
    inp = struct( 'analind', complvec(1), ...
                    'lfl', 1, ...
                    'extraL', [], ...
                    'ovrwrt', 2, ...
                    'multiflag', 2, ...
                    'saveparam', 2, ...
                    'saveCV1', 2, ...
                    'loadparam', 2, ...
                    'batchflag', 2);
end
na_str = '?'; inp.datatype = 'VISdatamat'; 
% Resolves bug when running in batch mode
if ~isfield(inp,'extraL') , inp.extraL=[]; end
OverWriteStr = []; GridSelectStr = []; LoadModelsStr = []; LoadParamsStr = []; LoadStr = []; SaveStr = []; MultiStr = []; ExtraLStr = []; SaveCV1Str = [];
OverWriteAct = []; GridSelectAct = []; LoadModelsAct = []; LoadParamsAct = []; LoadAct = []; SaveAct = []; MultiAct = []; ExtraLAct = []; SaveCV1Act = [];

%% Configure menu
if isfield(inp,'analind'), AnalSelStr = sprintf('Analysis %g', inp.analind); else, AnalSelStr = na_str;  end
AnalSelectStr = sprintf('Choose analysis to work on [ %s ]|', AnalSelStr);                                                  AnalSelectAct = 1;         

analysis = dat.analysis{inp.analind}; 
[inp.permfl, inp.permmode, inp.permsig] = get_permflag( analysis );

if ~isempty(analysis)
    
    % Initialize global parameters for the selected analysis
    nk_SetupGlobVars2(analysis.params, 'setup_main', 0); 
    MULTI = analysis.params.TrainParam.MULTI;
    
    % Compute from scratch or use pre-computed datamats ?
    LFL_opts        = {'Compute from scratch',sprintf('Use precomputed %s',inp.datatype)};                                      
    ModeStr         = sprintf('Operation mode of visualization module [ %s ]|',LFL_opts{inp.lfl});                          ModeAct = 2;
    
    if inp.lfl == 1,
        % from scratch
        OVRWRT_opts     = {'Overwrite existing','Do not overwrite'};       
        OverWriteStr = sprintf('Overwrite existing %s files [ %s ]|', inp.datatype, OVRWRT_opts{inp.ovrwrt}) ;              OverWriteAct = 3; 
    else
        % precomputed
        nVisFiles = na_str;
        if isfield(inp,'vismat') && ~isempty(inp.vismat), 
            selGrid = ~cellfun(@isempty,inp.vismat); inp.GridAct = selGrid;
            nVisFiles = sprintf('%g selected', sum(selGrid(:))); 
        end     
        OverWriteStr = sprintf('Specify %s files [ %s ]|', inp.datatype, nVisFiles);                                        OverWriteAct = 3; 
    end
    
    % Retrieve CV2 partitions to operate on
    if ~isfield(inp,'GridAct'), inp.GridAct = analysis.GDdims{1}.GridAct; end;                                              
    GridSelectStr = sprintf('Select CV2 partitions to operate on [ %g selected ]|',  sum(inp.GridAct(:)));                  GridSelectAct = 4;
        
    % Check multi-group settings
    if ~isempty(MULTI) && MULTI.flag && ~MULTI.BinBind  
        MULTI_opts      = {'Compute at multi-group optimum','Compute at binary optima'};                                        
        MultiStr = sprintf('Visualize patters in the multi-group setting [ %s ]|',MULTI_opts{inp.multiflag});               MultiAct = 5;
    end

    % Configure extra label dialogue
    if inp.permfl
        if ~isempty(inp.extraL)
            m = size(inp.extraL.L,2); 
            ExtraLStr = sprintf('Deactivate assessment of prognostic generalization [ %g extra label(s) defined ]|', m);     ExtraLAct = 11;
        else
            ExtraLStr = sprintf('Add extra labels for the assessment of prognostic generalization|');                        ExtraLAct = 10;
        end
    end 

    % Configure loading of pre-existing parameters and models
    if inp.saveparam == 2 && inp.lfl == 1
        LOAD_opts        = {'yes', 'no'}; 
        LoadStr = sprintf('Use saved pre-processing params and models [ %s ]|', LOAD_opts{inp.loadparam});                  LoadAct = 7;
        if inp.loadparam == 1
            if isfield(inp,'optpreprocmat'), 
                selGridPreproc = ~cellfun(@isempty,inp.optpreprocmat);
                nParamFiles = sprintf('%g files selected', sum(selGridPreproc(:))); 
            else, 
                nParamFiles = na_str; 
            end
            LoadParamsStr = sprintf('Select preprocessing parameter files [ %s ]|' ,nParamFiles);                           LoadParamsAct = 8;
            if isfield(inp,'optmodelmat'), 
                selGridModel = ~cellfun(@isempty,inp.optmodelmat);
                nModelFiles = sprintf('%g files selected', sum(selGridModel(:))); 
            else, 
                nModelFiles = na_str; 
            end
            LoadModelsStr = sprintf('Select model files [ %s ]|',nModelFiles);                                              LoadModelsAct = 9;
        end
    end
    
    % If loading of pre-existing models and params is not chosen, allow to
    % save the computed params and models to disk
    if inp.loadparam == 2 && inp.lfl == 1
        SAVE_opts       = {'yes', 'no'};   
        SaveStr = sprintf('Save pre-processing params and models to disk [ %s ]|', SAVE_opts{inp.saveparam});               SaveAct = 6;
        if inp.saveparam == 1
            SaveCV1Str = sprintf('Save pre-processing params at CV1 level [ %s ]|', SAVE_opts{inp.saveCV1});                SaveCV1Act = 12;
        end
    end

end

%% Build interactive menu
menustr = [ AnalSelectStr ...
           ModeStr ...
           ExtraLStr ...
           OverWriteStr ...
           GridSelectStr ...
           MultiStr ...
           SaveStr ...
           SaveCV1Str ...
           LoadStr ...
           LoadParamsStr ... 
           LoadModelsStr ];

menuact = [ AnalSelectAct ...
            ModeAct ...
            ExtraLAct ...
            OverWriteAct ...
            GridSelectAct ...
            MultiAct ...
            SaveAct ...
            SaveCV1Act ...
            LoadAct ...
            LoadParamsAct ...
            LoadModelsAct ];       

disallow = false;

%% Check whether all parameters are available
if ~sum(inp.GridAct(:)) || isempty(inp.analind), disallow = true; end

if inp.loadparam == 1
    if ~isfield(inp,'optpreprocmat') || isempty(inp.optpreprocmat), disallow = true; end
    if ~isfield(inp,'optmodelmat') || isempty(inp.optmodelmat), disallow = true; end
end

if ~disallow, menustr = [menustr '|PROCEED >>>']; menuact = [menuact 99]; end

%% Display menu and act on user selections
nk_PrintLogo
mestr = 'Visualization module run-time configuration'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr);
if inp.batchflag == 2, act = nk_input(mestr, 0, 'mq', menustr, menuact); end

switch act
    case 0
        return
    case 1
        showmodalvec = []; analind = inp.analind; 
        if length(dat.analysis)>1, t_act = 1; brief = 1;
            while t_act>0, 
                [t_act, analind, ~, showmodalvec, brief ] = nk_SelectAnalysis(dat, 0, navistr, analind, [], 1, showmodalvec, brief); 
            end
            if ~isempty(analind), inp.analind = analind ; end
        end
        inp.GridAct = dat.analysis{inp.analind}.GDdims{1}.GridAct;
        
    case 2
        lfl = nk_input('Define run-time mode of visualization module',0,'mq',strjoin(LFL_opts, '|'),[1,2],inp.lfl);
        if lfl, inp.lfl = lfl; end        
    case 3
        switch inp.lfl
            case 1
                if inp.ovrwrt == 1, inp.ovrwrt = 2; elseif inp.ovrwrt  == 2, inp.ovrwrt = 1; end
            case 2
                tdir = create_defpath(dat.analysis{inp.analind});
                inp.vismat = nk_GenDataMaster(dat.id, 'VISdatamat', CV, [], tdir);
        end
    case 4
        [operms,ofolds] = size(CV.TrainInd);
        tact = 1; while tact > 0 && tact < 10, [ tact, inp.GridAct ] = nk_CVGridSelector(operms, ofolds, inp.GridAct, 0); end
    case 5
        if inp.multiflag  == 1, inp.multiflag = 2; elseif inp.multiflag == 2,  inp.multiflag = 1; end
    case 6
        if inp.saveparam == 1, inp.saveparam = 2; elseif inp.saveparam == 2,  inp.saveparam = 1; end
    case 7
        if inp.loadparam == 1, inp.loadparam = 2; elseif inp.loadparam == 2,  inp.loadparam = 1; end
    case 8
        tdir = create_defpath(dat.analysis{inp.analind});
        optpreprocmat = nk_GenDataMaster(dat.id, 'OptPreprocParam', CV, [], tdir);
        if ~isempty(optpreprocmat), inp.optpreprocmat = optpreprocmat; end
    case 9
        tdir = create_defpath(dat.analysis{inp.analind});
        optmodelmat = nk_GenDataMaster(dat.id, 'OptModel', CV, [], tdir);
        if ~isempty(optmodelmat), inp.optmodelmat = optmodelmat; end
    case 10
        Ldef = []; Lsize = Inf; LNameDef = []; ActStr = 'Define'; 
        if isfield(inp,'extraL')
            if isfield(inp.extraL,'L')
                Ldef = inp.extraL.L; 
                Lsize = size(LDef,2); 
                ActStr = 'Modify'; 
            end
        end
        inp.extraL.L = nk_input([ ActStr ' extra label matrix for permutation testing'],0,'e',Ldef,[numel(NM.label),Lsize]);
        if isfield(inp.extraL,'Lnames'), LNameDef = inp.extraL.Lnames; end
        inp.extraL.Lnames = nk_input([ ActStr ' cell array of string descriptors for extra labels'],0,'e',LNameDef,[1 Lsize]);
    case 11
        del_extraL = nk_input('Do you really want to delete the extra labels?',0,'yes|no',[1,0]);
        if del_extraL, inp.extraL = []; end
    case 12
        if inp.saveCV1 == 1, inp.saveCV1 = 2; elseif inp.saveCV1 == 2,  inp.saveCV1 = 1; end
    case 99
        NM.runtime.curanal = inp.analind;
        dat.analysis{inp.analind} = VisModelsPrep(dat, inp, analysis);
end

function tdir = create_defpath(analysis)
 
rootdir = analysis.rootdir;
algostr = getAlgoStr(analysis);
VIS = analysis.params.TrainParam.VIS{1};
        
if VIS.PERM.flag
    switch VIS.PERM.mode 
        case 1
            visdir = 'VISUAL_L';
        case 2
            visdir = 'VISUAL_F';
        case 3
            visdir = 'VISUAL_LF';
    end
else
    visdir = 'VISUAL';
end
tdir = fullfile(rootdir, algostr, visdir);


function analysis = VisModelsPrep(dat, inp1, analysis)
global CV COVAR MODEFL FUSION MULTILABEL

if inp1.multiflag   == 2, inp1.multiflag    = 0; end
if inp1.saveparam   == 2, inp1.saveparam    = 0; end
if inp1.loadparam   == 2, inp1.loadparam    = 0; end
if inp1.ovrwrt      == 2, inp1.ovrwrt       = 0; end

F = 1; nF = 1;
if ~isempty(FUSION)        
    F = analysis.params.TrainParam.FUSION.M;
    nF = numel(F); if FUSION.flag < 3, nF = 1; end
end

switch inp1.lfl
    case 1
        % *********** CONSTRUCT VISUALIZATION ANALYSIS INPUT STRUCTURE ************
        inp1.analmode    = 0;
        inp1.covstr = '';
        if ~isempty(COVAR)
            inp1.covstr = dat.covnames{1};
            for j=2:length(COVAR)
                inp1.covstr = [inp1.covstr ', ' dat.covnames{j}];
            end
        end        
    case 2
        inp1.analmode = 1;
        inp1.covstr = '';
end

inp1.nclass = 1;if strcmp(MODEFL,'classification'), inp1.nclass = numel(CV.class{1,1}); end
if inp1.permfl
    permmodestr = {'L','F','LF'};
    permstr = sprintf('VISUAL_%s',permmodestr{inp1.permmode});
else
    permstr = 'VISUAL';
end
if isfield(analysis,'rootdir') && exist(analysis.rootdir,'dir')
    inp1.rootdir = fullfile(analysis.rootdir,analysis.params.TrainParam.SVM.prog,permstr);
else
    inp1.rootdir = fullfile(pwd,analysis.params.TrainParam.SVM.prog,permstr);
end
inp1.procdir = fullfile( analysis.rootdir, 'proc');
if ~exist(inp1.rootdir,'dir'), mkdir(inp1.rootdir); end

%%%%%%%%%%%%%%%%%%%%%%% RUN VISUALIZATION ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nF
    inp2            = nk_SetFusionMode2(dat, analysis, F, nF, i);
    inp             = catstruct(inp1,inp2); clear inp2;
    inp.curlabel    = 1;
    VIS = analysis.params.TrainParam.VIS{analysis.params.TrainParam.FUSION.M(i)};
    % Per default do not activate normalization of the weight vector
    % (30.12.2018)
    if isfield(VIS,'norm'), inp.norm = VIS.norm; else, inp.norm = false; end

    for j = 1:MULTILABEL.dim
        if MULTILABEL.flag && MULTILABEL.dim>1, 
            fprintf('\n\n');cprintf('*black','====== Working on label #%g ====== ',j); inp.curlabel = j; 
        end
        analysis.visdata{i,j} = nk_VisModels(inp, dat.id, inp.GridAct);
    end
end

function [ permfl, permmode, permsig] = get_permflag ( analysis )

permfl = false; permmode =[]; permsig = [];
varind = analysis.params.TrainParam.FUSION.M; 
for i=1:numel(varind)
    if isfield(analysis.params.TrainParam.VIS{varind(i)},'PERM') && analysis.params.TrainParam.VIS{varind(i)}.PERM.flag
        permfl = true; 
        permmode = analysis.params.TrainParam.VIS{varind(i)}.PERM.mode;
        permsig = analysis.params.TrainParam.VIS{varind(i)}.PERM.sigflag;
        break;
    end
end

