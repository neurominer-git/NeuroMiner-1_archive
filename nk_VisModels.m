function visdata = nk_VisModels(inp, id, GridAct, batchflag)
% =========================================================================
% visdata = nk_VisModels(inp, strout, id, GridAct, batchflag)
% =========================================================================
% Visualisation module of NeuroMiner
% ...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2017

global SVM SAV RFE MODEFL CV VERBOSE FUSION MULTILABEL EVALFUNC

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
visdata         = [];                               % Initialize with empty output
switch inp.analmode
    case 0
        ovrwrt  = inp.ovrwrt;                       % overwrite existing data files
    case 1
        vismat  = inp.vismat;                       % Visualization datamat
end
saveparam       = inp.saveparam;
loadparam       = inp.loadparam;
strout          = inp.varstr;                       % Suffix for filename indicating modality assignment of file
nclass          = inp.nclass;                       % Number of binary comparisons
analysis        = inp.analysis;                     % GDanalysis structure to be used
varind          = inp.tF;                           % Current modality index
ol              = 0;                                % Counter for computing mean vectors
ll              = 1;                                % Counter for looping through CV2 arrays
[ix, jx]        = size(CV.TrainInd);                % No. of Perms / Folds at CV2 level
ind0            = false(1,ix*jx);
algostr         = GetMLType(SVM);                   % Algorithm descriptions
FullPartFlag    = RFE.ClassRetrain;                 % Has the user activated full CV1 partition retraining ?
nM              = numel(inp.tF);                    % Number of modalities with independent preprocessing
decompfl        = false(1,nM);                      % flag for factorization methods during preprocessing
permfl          = false(1,nM);                      % flag for permutation mode 
sigfl           = false(1,nM);                      % flag for significance estimation mode
nperms          = ones(1,nM);                       % number of permuations per modality
pmode           = ones(1,nM);                       % Permutation mode

% Loop through modalities (in early fusion: nM = 1)
Dall = 0;

% Check whether you have to run label imputation
IMPUTE.flag = false;
if iscell(inp.PREPROC), iPREPROC = inp.PREPROC{1}; else iPREPROC = inp.PREPROC; end    
if isfield(iPREPROC,'LABELMOD') && isfield(iPREPROC.LABELMOD,'LABELIMPUTE'); 
    IMPUTE = iPREPROC.LABELMOD.LABELIMPUTE; 
    IMPUTE.flag = true; 
end
linsvmfl = determine_linsvm_flag(SVM);

% Always set to binary preprocessing (unless true multi-group learners have
% been intergrated in NM)
BINMOD = iPREPROC.BINMOD; 
clc

inp.id = id;

for i = 1 : nM
    
    % Dimensionality of current modality
    D = getD(FUSION.flag, inp, i);
    
    % Dimensionality of the (concatenated feature space)
    Dall = Dall + D;
    
    % Activate preprocessing params of current modality
    switch FUSION.flag
        case 2
            iPREPROC = inp.PREPROC{i}; 
            iVis = inp.VIS{i};
        otherwise
            iPREPROC = inp.PREPROC;
            iVis = inp.VIS;
    end
    
    % Determine if factorization methods are involved in current
    % preprocessing chain
    if isfield(iPREPROC,'ACTPARAM')
        for zu = numel(iPREPROC.ACTPARAM) : -1 : 1
            if strcmp(iPREPROC.ACTPARAM{zu}.cmd,'reducedim') ...
                || strcmp(iPREPROC.ACTPARAM{zu}.cmd,'remvarcomp') 
                decompfl(i) = true; 
            end
        end
    end
    
    % Check whether permutation mode is activated in the visualisation setup of current modality 
    if isfield(iVis,'PERM') && iVis.PERM.flag
        if ~isfield(iVis.PERM,'mode')
            pmode(i) = 1; 
        else
            pmode(i) = iVis.PERM.mode;
        end
        permfl(i)      = true; 
        nperms(i)      = iVis.PERM.nperms;
        compfun        = nk_ReturnEvalOperator(SVM.GridParam);
        
        % Generate or load permutation matrix from file.
        cprintf('black*','\nPERMUTATION MODE ENABLED.')
        if  pmode(i)==1
            permfile = fullfile(inp.rootdir,[SAV.matname '_VISpermmat_ID' id '.mat']);
            if ~exist(permfile,'file')
                fprintf('\nCreating parent permutation matrix with %g perms', nperms(1))
    %             if nclass>1
    %                 indpermA = cell(nclass,1);
    %                 for curclass=1:nclass
    %                     fprintf('\n\tBinary comparison: %g', curclass);
    %                     groups = CV.class{1,1}{curclass}.groups;
    %                     if numel(groups)==2
    %                         N = sum(inp.labels == groups(1) | inp.labels == groups(2)) + sum(isnan(inp.labels));
    %                     else
    %                         N = size(inp.labels,1);
    %                     end
    %                     indpermA{curclass} = nk_VisXPermHelper('genpermlabel', N, nperms(1));
    %                 end
    %             else
                    indpermA = nk_VisXPermHelper('genpermlabel', size(inp.labels,1), nperms(1));
     %           end
                save(permfile,'indpermA');
            else
                fprintf('\nLoading parent permutation matrix from file: \n%s', permfile);
                load(permfile,'indpermA');
                if size(indpermA,2)~= nperms(1)
                    error('Number of permutations in matrix (%g) does not match number of expected permutations (%g)', size(indpermA,2), nperms(1))
                end
                if size(indpermA,1) ~= size(inp.labels,1) 
                    error('Number of cases in permutation matrix (%g) does not match number of cases in NM file (%g)', size(indpermA,1), size(inp.labels,1))
                end
            end
        end
        % Check whether reconstruction of only those components that are significant is
        % activated by the user.
        if decompfl(i) && isfield(iVis.PERM,'sigflag') && iVis.PERM.sigflag,
            sigfl(i) = true;
            sigflFDR = false;
        end
    end
    if isfield(iVis,'use_template') && iVis.use_template,
        fprintf('\nPerform template processing!')
        templateflag(i) = true;
    end        
end

% Initialize CV2 data containers
I = nk_VisXHelper('init', nM, nclass,  decompfl, permfl, ix, jx);
if any(permfl), 
    I.VCV2MPERM_S = cell(size(inp.labels,1),nclass,nperms(1)); 
    I.VCV2MORIG_S = cell(size(inp.labels,1),nclass); 
end

% Obtain feature labels for the selected modalities
featnames = []; if isfield(inp,'featnames'), featnames = inp.featnames; end

% Multi-Group processing ?
multiflag = false; if isfield(inp,'multiflag') && ~isempty(inp.multiflag), multiflag = inp.multiflag; end
if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end
multlabelstr = '';  if MULTILABEL.flag, multlabelstr = sprintf('t%g',inp.curlabel); end

 % Do we have to scale the labels?
[ inp ] = nk_ApplyLabelTransform( inp.PREPROC, MODEFL, inp );

% %%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds
        
        fprintf('\n---------------------------------------------------------------------------------------------')
        if ~GridAct(f,d), 
            ll=ll+1;
            fprintf('\nSkipping CV2 [%g,%g] (user-defined).',f,d)
            continue 
        end;
        
        [iy, jy] = size(CV.cvin{f,d}.TrainInd); % No. of Perms / Folds at CV1 level
        
        operm = f; ofold = d;
        oVISpath = nk_GenerateNMFilePath(inp.rootdir, SAV.matname, inp.datatype, multlabelstr, strout, id, operm, ofold);
        OptModelPath = nk_GenerateNMFilePath( inp.rootdir, SAV.matname, 'OptModel', [], inp.varstr, inp.id, operm, ofold);
    
        switch inp.analmode 
            
            case 0
        
                %%%%%%%%%%%%%%%%%%%%%%%%% USE PRECOMPUTED ? %%%%%%%%%%%%%%%%%%%%%%%
                
                if exist(oVISpath,'file') && ~ovrwrt && ~batchflag
                    
                    [~, onam] = fileparts(oVISpath);
                    fprintf('\nVISdatamat found for CV2 [%g,%g]:',f,d)
                    fprintf('\nLoading: %s',onam)
                    [I, I1] = nk_VisXHelper('accum', nM, nclass, decompfl, permfl, ix, jx, I, inp, ll, nperms, oVISpath);
                    if isfield(I1,'PCV1SUM'), PCV1SUMflag=true; else PCV1SUMflag = false; end
                    if any(permfl)
                        for h=1:nclass
                            switch MODEFL
                                case 'classification'
                                    TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                                case 'regression'
                                    TsInd2 = CV.TestInd{f,d};
                            end
                            I.VCV2MORIG_S(TsInd2,h) = cellmat_mergecols(I.VCV2MORIG_S(TsInd2,h), num2cell(I1.DS{h},2));
                            for perms = 1:nperms(1), I.VCV2MPERM_S(TsInd2,h,perms) = cellmat_mergecols(I.VCV2MPERM_S(TsInd2,h,perms), num2cell(I1.DS_perm{h}(:,:,perms),2)); end
                        end
                    end
                    ll=ll+1;
                    ind0(ll) = true;
                    ol=ol+1; continue
                    
                elseif exist(oVISpath,'file') && batchflag
                    
                    % in batch mode we do not compute statistics across the
                    % CV2 partitions
                    [~, onam] = fileparts(oVISpath);
                    fprintf('\nVISdatamat found for CV2 [%g,%g]:\n%s',f,d,onam)
                    fprintf('\nBatch mode detected. Continue.')
                    ll=ll+1;
                    continue

                end

                %%%%%%%%% GET PREPROCESSING PARAMETERS FOR CUR. CV2 PART. %%%%%%%%%
                % First generate parameter array for preprocessing based on
                % the trained base learners in the ensemble. This saves
                % computational resources because we are not going to preprocess the 
                % data with all possible parameter combinations specified by the user, 
                % but only with those chosen by the NM training process.
                
                inp.f = f; inp.d = d; inp.ll = ll;  
                % Parameter flag structure for preprocessing
                paramfl = struct('use_exist', loadparam, ...
                                 'found', false, ...
                                 'write', 1, ...
                                 'multiflag', multiflag);
                
                % Apply prerpocessing on the entire data and use these
                % parameters to adjust for arbitrary PCA rotations through 
                % the Procrustes transform 
                if exist('templateflag','var') && any(templateflag), paramfl.templateflag = true; end
                               
                % Compute params
                [ contfl, analysis, mapY, GD, ~, Param, paramfl ] = nk_ApplyTrainedPreproc2(analysis, inp, paramfl);
                if contfl, continue; end
                    
                % Prepare containers & initialize matrices for CV1-level
                % weight vector relevance metrics
                % =========================================================
                ol                                = ol+1;
                [~, I1] = nk_VisXHelper('initI1', nM, nclass, decompfl, permfl);
                GDFEAT                         = GD.FEAT; 
                GDVI                           = GD.VI; 
                if inp.stacking
                    mChnl = GD.nM_cnt;
                end
                clear GD
                
                % Try to load models from disk if user
                % chose to do this
                fndMD = false; 
                if loadparam && isfield(inp,'optmodelmat') && exist(inp.optmodelmat{operm,ofold},'file')
                    fprintf('\nLoading OptModel: %s', inp.optmodelmat{operm,ofold});
                    load(inp.optmodelmat{operm,ofold},'MD'); fndMD = true; 
                end
                if ~fndMD, MD = cell(nclass,1); end
                
                % ---------------------------------------------------------
                if ~VERBOSE,fprintf('\n\nComputing visualizations'), end
                
                if permfl, 
                    I1.TS           = cell(nclass,1);
                    I1.DS           = cell(nclass,1);
                    I1.TS_perm      = cell(nclass,1);
                    I1.DS_perm      = cell(nclass,1);
                end
                
                for h=1:nclass % Loop through binary comparisons
                    
                    if nclass > 1, 
                        fprintf('\n*** %s #%g ***',algostr, h);                            
                    end
                    
                    TsInd = mapY.TsInd{h}; 
                    switch MODEFL
                        case 'classification'
                            TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                        case 'regression'
                            TsInd2 = CV.TestInd{f,d};
                    end
                    %% Step 1: Get optimal model parameters
                    % Retrieve optimal parameters from precomputed analysis structure
                    % Differentiate according to binary or multi-group mode
                    [Ps, Pspos, nP, Pdesc] = nk_GetModelParams2(analysis, multiflag, ll, h);
                    
                    % Computation of feature selection probabilities works
                    % only if features used for classification live in the
                    % original feature space and not in decomposed spaces
                    % (e.g. PCA, NMF-derived features)
                    
                    % Allocate memory to store CV1 ensemble patterns
                    il = 1; kil=1; ill = getModelNumDim(h,iy,jy,nP,Pspos,GDFEAT);
                    fprintf('\nNeed to evaluate %g models in this CV2 partition', ill);
                    
                    if permfl, 
                        I1.VCV1WPERM{h} = cell(ill,1); 
                        I1.VCV1MPERM{h} = nan(ill,1); 
                    end
                    
                    for n=1:nM
                        
                        % Retrieve dimensionality of target space
                        D                          = getD(FUSION.flag, inp, n);
                        % Setup container for weight storage
                        I1.VCV1{h,n}               = zeros(D, ill,'single'); 
                        I1.VCV1SUM{h,n}            = zeros(D, 1,'single');
                        I1.VCV1SQ{h,n}             = I1.VCV1SUM{n};
                        I1.VCV1MEAN{h,n}           = I1.VCV1SUM{n};
                        I1.VCV1STD{h,n}            = I1.VCV1SUM{n};
                        
                        % Prepare for analysis without factorization
                        if any(~decompfl), 
                            I1.PCV1SUM{h, n}                    = zeros(D, 1,'single'); 
                            I1.VCV1PEARSON{h, n}                = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1SPEARMAN{h, n}               = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1PEARSON_UNCORR_PVAL{h, n}    = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}   = nan(D, iy*jy*nP,'single');
                            I1.VCV1PEARSON_FDR_PVAL{h, n}       = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1SPEARMAN_FDR_PVAL{h, n}      = nan(D, iy*jy*nP,'single'); 
                        end
                        
                        % Prepare for permutation analysis
                        if permfl
                            nTs                 = size(TsInd,1);
                            I1.TS{h}            = nan(nTs, ill);
                            I1.DS{h}            = nan(nTs, ill);
                            I1.DS_perm{h}       = nan(nTs, ill, nperms(1));
                            I1.TS_perm{h}       = nan(nTs, ill, nperms(1));
                            I1.VCV1PERM{h, n}   = nan(D, ill,'single');
                            I1.VCV1PERM_FDR{h, n} = nan(D, ill,'single');
                            I1.VCV1ZSCORE{h, n} = nan(D, ill,'single');
                        end
                    end
                    
                    if ~fndMD , MD{h} = cell(nP,1); end
                    
                    for m = 1 : nP
                        
                        if nP>1, fprintf('\nExtracing model parameters at parameter node %g of %g', m, nP); end
                        % Prepare learning params
                        cPs = Ps(m,:); sPs = nk_PrepMLParams(Ps, Pdesc, m);
                        
                        % -----------------------------------------------------
                        % Construct pattern for every base learnern in
                        % current CV1 [k,l] partition:
                        %% CV1 LOOP
                        P_str = nk_DefineMLParamStr(cPs, analysis.Model.ParamDesc, h);
                        
                        if ~fndMD,MD{h}{m} = cell(iy,jy); end
                        
                        for k=1:iy % permutations

                            for l=1:jy % folds
                                
                                Fkl = GDFEAT{Pspos(m)}{k,l,h}; 
                                    
                                % Determine number of features in mask and
                                % convert feature mask to logical index
                                ul=size(Fkl,2);
                                if ~islogical(Fkl),F = Fkl ~= 0; else F = Fkl; end

                                if VERBOSE
                                    fprintf('\n');cprintf('*black',['Constructing predictive pattern(s) in CV2 [%2g ,%2g], ' ...
                                    'CV1 [%2g ,%2g]: %g model(s), %s ML params [ %s ]. '], f, d, k, l, ul, algostr, P_str); 
                                else
                                    fprintf('\n');cprintf('*black','CV2 [%2g, %2g ], CV1 [ %2g, %2g ]: %g model(s) ',f, d, k, l, ul) ;
                                end

                                CVInd   = mapY.CVInd{k,l}{h};
                                TrInd   = mapY.TrInd{k,l}{h};
                                                        
                                % Set the pointer to the correct mapY shelf
                                for n=1:numel(paramfl)
                                    pnt = 1;
                                    if ~BINMOD
                                         if isfield(paramfl{n},'PREPROC') && ...
                                           isfield(paramfl{n},'PXfull') && ...
                                           ~isempty(paramfl{n}.P{1})
                                            pnt = m;
                                            break   
                                        end
                                    else
                                        if isfield(paramfl{n},'PREPROC') && ...
                                           isfield(paramfl{n},'PXfull') && ...
                                           ~isempty(paramfl{n}.P{h})
                                            pnt = m;
                                            break   
                                        end
                                    end
                                end
                                
                                %%%%%% RECOMPUTE ORIGINAL MODEL %%%%%%
                                % get CV1 training and test data
                                if BINMOD, hix = h; else, hix =1; end
                                [ modelTr , modelCV, modelTs] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, [], Param{1}(k,l,hix), pnt); 
                                switch FUSION.flag
                                    case 2
                                        for n=1:nM, 
                                            [~,~,~,~, ParamX{n} ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, [], Param{n}(k,l,hix), pnt); 
                                        end
                                    otherwise
                                        ParamX = Param{1}(k,l,hix).TrainedParam;
                                end      

                                % Get label info
                                modelTrL = mapY.TrL{k,l}{h};
                                modelCVL = mapY.CVL{k,l}{h};
                                
                                % Impute labels if needed
                                [modelTrL, modelTr, TrInd] = nk_LabelImputer( modelTrL, modelTr, TrInd, sPs, IMPUTE);
                                
                                % Concatenate Training and CV data if needed
                                if FullPartFlag, 
                                    modelTr = [modelTr; modelCV ]; 
                                    modelTrL = [modelTrL; modelCVL]; 
                                    TrInd = [TrInd; CVInd]; 
                                end
                                modelTr = modelTr(TrInd,:);
                                
                                % Prepare permutation operations
                                if any(permfl),
                                    indperm = []; indpermfeat = [];
                                    switch MODEFL
                                        case 'classification'
                                            pTrInd = CV.TrainInd{f,d}(CV.class{f,d}{h}.TrainInd{k,l}); 
                                            pCVInd = CV.TrainInd{f,d}(CV.class{f,d}{h}.TestInd{k,l}); 
                                        otherwise
                                            pTrInd = CV.TrainInd{f,d}(CV.cvin{f,d}.TrainInd{k,l}); 
                                            pCVInd = CV.TrainInd{f,d}(CV.cvin{f,d}.TestInd{k,l});
                                    end
                                    if size(pTrInd,2)>1,pTrInd=pTrInd'; end
                                    if size(pCVInd,2)>1,pCVInd=pCVInd'; end
                                    if FullPartFlag, pTrInd = [pTrInd; pCVInd];end
                                    if pmode(1)==1, indperm = indpermA(pTrInd,:); end
                                    modelTs = modelTs(TsInd,:); 
                                    modelTsL = mapY.TsL{h};
                                else
                                    modelTs = []; modelTsL = [];
                                end
                                
                                if ~fndMD, MD{h}{m}{k,l} = cell(ul,1); end
                                
                                % Loop through base learners' feature masks
                                for u=1:ul
                                     
                                    % Extract features according to mask
                                    try
                                        Ymodel = nk_ExtractFeatures(modelTr, F, [], u);
                                    catch
                                        error('Dimensionality mismatch between training data matrix and feature selection mask. Check your settings')
                                    end
                                    Find = F(:,u);
                                    %If permutation mode expects feature
                                    % permutation prepare for this here:
                                    if any(permfl) && pmode(1) > 1, indpermfeat = nk_VisXPermHelper('genpermfeats', sum(F(:,u)), nperms(1), modelTrL); end
                                    
                                    if ~isempty(GDVI{Pspos(m)})
                                        Vind = GDVI{Pspos(m)}{k,l,h}(:,u);
                                    else
                                        Vind = true(size(Find,1),1);
                                    end
                                    
                                    % Model computation
                                    if ~fndMD, [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(Ymodel, modelTrL, sPs, 1); end
                                    
                                    if inp.stacking
                                        vec_mj = [];
                                        for mj = 1:inp.nD
                                            vec_mj = [vec_mj; mj*ones(mChnl(mj),1)];
                                        end
                                    end
                                    
                                    if ~any(permfl)
                                        % Compute original weight map in input space
                                        [Tx, Psel, Rx, SRx, ~, PAx] = nk_VisXWeight2(inp, MD{h}{m}{k,l}{u}, Ymodel, modelTrL, varind, ParamX, Find, Vind, decompfl, pnt);
                                    else
                                        fprintf('\n\t%3g | OptModel =>', u);
                                         % Build permutation index array if
                                        % needed and compute further variables
                                        % needed for perm-based stats
                                        [perf_orig, I1.TS{h}(:,il), I1.DS{h}(:,il)] = nk_GetTestPerf(modelTs, modelTsL, Find, MD{h}{m}{k,l}{u}, modelTr); 
                                        fprintf(' %1.2f',perf_orig)
                                        % if sigfl = true
                                        % Determine significant components
                                        % through non-parametric permutation
                                        if sigfl
                                            fprintf(' | Significant pattern components (%g perms):\t',nperms(1));
                                            % Original weight vector
                                            [~,~,~,~, Vx] = nk_VisXWeight2(inp, MD{h}{m}{k,l}{u}, Ymodel, modelTrL, varind, ParamX, Find, Vind, decompfl, pnt, false);
                                            Vx_perm = zeros( size(Vx,1), nperms(1) );
                                            MD_perm = cell(nperms(1),1);
                                            for perms = 1:nperms(1)
                                                fprintf('+');
                                                % Train permuted model
                                                [L_perm, Ymodel_perm]       = nk_VisXPermY(Ymodel, inp.labels, pTrInd, pmode(1), indperm, indpermfeat, perms);
                                                [~, MD_perm{perms}]         = nk_GetParam2(Ymodel_perm, L_perm, sPs, 1);
                                                % Compute permuted weight vector in feature space
                                                [~,~,~,~,Vx_perm(:, perms)] = nk_VisXWeight2(inp, MD_perm{perms}, Ymodel_perm, L_perm, varind, ParamX, Find, Vind, decompfl, pnt, false);
                                            end
                                            % Determine significance by comparing observed weight vector
                                            % against null distribution of permuted models' weight vectors
                                            %(abs(w/(norm(w,2))
                                            I1.VCV1WPERM{il,h} = (sum(bsxfun(@ge,abs(Vx_perm),abs(Vx)),2)/nperms(1))/ul;
                                            if sigflFDR
                                                Fadd   = fdr_bh(I1.VCV1WPERM{il,h},0.5,'dep');
                                                FDRstr = '(FDR) ';
                                            else
                                                Fadd   = (I1.VCV1WPERM{il,h} <= 0.5);
                                                FDRstr = '(uncorr) ';
                                            end
                                            if ~sum(Fadd)
                                                [minP, Fadd] = min(I1.VCV1WPERM{il,h});
                                                fprintf('\tNo component significant at alpha %s = 0.5 => relaxing to max P = %g\n\t\t\t\t\t\t\t\t\t\t',FDRstr, minP );
                                            else
                                                fprintf('\t%g / %g components significant at alpha %s = 0.5\n\t\t\t\t\t\t\t\t\t\t',sum(Fadd), numel(Fadd),FDRstr);
                                            end
                                        else
                                            Fadd = true(size(F,1),1);
                                        end
                                        
                                        fprintf(' | Permuting:\t');
                                        % Compute original weight map in input space
                                        [Tx, Psel, Rx, SRx, ~, PAx ] = nk_VisXWeight2(inp, MD{h}{m}{k,l}{u}, Ymodel, modelTrL, varind, ParamX, Find, Vind, decompfl, pnt, [], Fadd);
                                        Tx_perm = cell(1,nM); Px_perm = zeros(1,nperms(1));
                                        for n=1:nM, Tx_perm{n} = zeros(size(Tx{n},1),nperms(n)); end
                                        for perms = 1:nperms(1)
                                            if ~sigfl, 
                                                % Train permuted model
                                                [ L_perm, Ymodel_perm ] = nk_VisXPermY(Ymodel, inp.labels, pTrInd, pmode(1), indperm, indpermfeat, perms);
                                                [~, MDs] = nk_GetParam2(Ymodel_perm, L_perm, sPs, 1);
                                            else
                                                % Retrieve trained permuted model
                                                MDs = MD_perm{perms};
                                            end
                                            % Compute permuted model test performance
                                            [perf_perm, I1.TS_perm{h}(:,il,perms), I1.DS_perm{h}(:,il,perms)] = nk_GetTestPerf(modelTs, modelTsL, Find, MDs, modelTr);
                                            % Compare against original model performance
                                            if feval(compfun, perf_perm, perf_orig),
                                                fprintf('.'); 
                                                Px_perm(perms) = Px_perm(perms) + 1;
                                            end
                                            % Compute permuted weight map in input space
                                            TXperms = nk_VisXWeight2(inp, MDs, Ymodel_perm, L_perm, varind, ParamX, Find, Vind, decompfl, pnt, [], Fadd);
                                            for n=1:nM, Tx_perm{n}(:,perms) = TXperms{n}; end
                                        end
                                        
                                        % Model significance
                                        I1.VCV1MPERM{h}(il) = (sum(Px_perm) / nperms(1));
                                        
                                        fprintf('\t');cprintf('*black', ' P=%1.3f ', I1.VCV1MPERM{h}(il));
                                        
                                        % Loop through modalities
                                        for n=1:nM
                                            % Define index vector to
                                            % original space of modality
                                            Fpind = any(Tx_perm{n},2)'; 
                                            
                                            % Now compute the P value vector:
                                            Pvals = sum(bsxfun(@ge,abs(Tx_perm{n}(Fpind,:)),abs(Tx{n}(Fpind))),2)/nperms(1);
                                            
                                            % ... and the Z score vector:
                                            Zvals = bsxfun(@rdivide, Tx{n}(Fpind) - nanmean(Tx_perm{n}(Fpind,2)), nanstd(Tx_perm{n}(Fpind,:),[],2));
                                            
                                            if inp.stacking
                                                for mj = 1:inp.nD
                                                    mjPvals = mean(Pvals(vec_mj(Fpind) == mj));
                                                    mjC = mean(Zvals(vec_mj(Fpind) == mj));
                                                    I1.VCV1PERM{h,n}(mj,il) = mjPvals;
                                                    [~,~,~,I1.VCV1PERM_FDR{h,n}(mj, il)] = fdr_bh(mjPvals,0.05,'pdep');
                                                    I1.VCV1ZSCORE{h,n}(mj,il) = mjC;
                                                end
                                            else
                                                [ ~, ~, ~, badcoords ] = getD(FUSION.flag, inp, n); badcoords = ~badcoords;
                                                I1.VCV1PERM{h,n}(badcoords & Fpind,il) = Pvals;
                                                [~,~,~,I1.VCV1PERM_FDR{h,n}(badcoords & Fpind, il)] = fdr_bh(Pvals,0.05,'pdep'); 
                                                I1.VCV1ZSCORE{h,n}(badcoords & Fpind,il) = Zvals;
                                            end
                                            % and show how many uncorrected P
                                            % values are below alpha=0.05
                                            sigcomp = sum(I1.VCV1PERM{h,n}(:,il)<=0.05);
                                            sigmin  = min(I1.VCV1PERM{h,n}(:,il));
                                            FDRsigmin = min(I1.VCV1PERM_FDR{h,n}(:,il));
                                            fprintf('; Modality #%g: [ %g features <= 0.05, min P value (uncorr, FDR) = %1.6f, %1.6f ] ', n, sigcomp, sigmin, FDRsigmin);
                                        end
                                    end
                                    
                                    % Some additional computation if
                                    % factorization methods have not been used
                                    % in the preprocessing chain of a modality
                                    for n = 1:nM
                                        Fpind = any(Tx{n},2); 
                                        if inp.stacking
                                            for mj = 1:inp.nD
                                                Imj = vec_mj == mj & Fpind;
                                                I1.VCV1{h,n}(mj,il) = mean(Tx{n}(Imj));
                                                if ~decompfl(n) && u==1,  
                                                    I1.VCV1PEARSON{h, n}(mj,kil) = nanmean(Rx{n}(Imj));
                                                    I1.VCV1SPEARMAN{h, n}(mj,kil) = nanmean(SRx{n}(Imj));
                                                    I1.VCV1PEARSON_UNCORR_PVAL{h, n}(mj,kil) = nanmean(nk_PTfromR(Rx{n}(Imj), size(Ymodel,1), 2));
                                                    I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}(mj,kil) = nanmean(nk_PTfromR(SRx{n}(Imj), size(Ymodel,1), 2));
                                                    if linsvmfl
                                                        I1.VCV1PVAL_ANALYTICAL{h, n}(mj,kil) = nanmean(PAx{n}(Imj));
                                                    end
                                                end
                                            end
                                        else
                                            % Define index vector to
                                            % original space of modality
                                            [ ~, ~, ~, badcoords] = getD(FUSION.flag, inp, n); badcoords = ~badcoords;

                                            % Store results in CV1 container variables                                    
                                            % I1.numCV1parts(h, n) = I1.numCV1parts(h, n) + 1;
                                            I1.VCV1{h,n}(badcoords,il) = Tx{n};

                                            if ~decompfl(n) 
                                                if u==1,  
                                                    %% Compute univariate correlation coefficient for each feature
                                                    I1.VCV1PEARSON{h, n}(badcoords,kil) = Rx{n};
                                                    I1.VCV1SPEARMAN{h, n}(badcoords,kil) = SRx{n};
                                                    I1.VCV1PEARSON_UNCORR_PVAL{h, n}(badcoords,kil) = nk_PTfromR(Rx{n}, size(Ymodel,1), 2);
                                                    I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}(badcoords,kil) = nk_PTfromR(SRx{n}, size(Ymodel,1), 2);
                                                    [~,~,~,I1.VCV1PEARSON_FDR_PVAL{h, n}(badcoords,kil)] = fdr_bh(I1.VCV1PEARSON_UNCORR_PVAL{h, n}(badcoords,kil), 0.05, 'pdep');
                                                    [~,~,~,I1.VCV1SPEARMAN_FDR_PVAL{h, n}(badcoords,kil)] = fdr_bh(I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}(badcoords,kil), 0.05, 'pdep');
                                                end
                                                %% Comnpute feature selection probabilities %%
                                                if isempty(I1.PCV1SUM{h, n}), I1.PCV1SUM{h, n} = zeros(size(badcoords,2),1); end
                                                I1.PCV1SUM{h, n}(badcoords) = I1.PCV1SUM{h, n}(badcoords) + Psel{n}; 
                                                if linsvmfl
                                                     I1.VCV1PVAL_ANALYTICAL{h, n}(badcoords,kil) = PAx{n};
                                                     [~,~,~,I1.VCV1PVAL_ANALYTICAL_FDR{h, n}(badcoords,kil)] = fdr_bh(PAx{n},0.05,'pdep');
                                                end
                                            end
                                       end
                                    end
                                    il=il+1;
                                end
                                kil=kil+1;
                                clear Tx Vx Tx_perm Vx_perm tSRx tRx Rx SRx MD_perm
                                %fprintf(' Done.')
                            end
                        end  
                    end
                    if any(permfl)
                        % Compute CV2-level model significance
                        I1.VCV1MORIG_EVALFUNC_CV2{h} = feval(EVALFUNC, modelTsL, nm_nanmedian(I1.DS{h},2)); 
                        I1.VCV1MPERM_CV2{h} = zeros(nperms(1),1); 
                        I1.VCV1MPERM_EVALFUNC_CV2{h} = zeros(nperms(1),1);
                        I.VCV2MORIG_S(TsInd2,h) = cellmat_mergecols(I.VCV2MORIG_S(TsInd2,h), num2cell(I1.DS{h},2));
                        for perms = 1:nperms(1)
                            I1.VCV1MPERM_EVALFUNC_CV2{h}(perms) = feval(EVALFUNC, modelTsL, nm_nanmedian(I1.DS_perm{h}(:,:,perms),2)); 
                            I1.VCV1MPERM_CV2{h}(perms)          = feval(compfun, I1.VCV1MPERM_EVALFUNC_CV2{h}(perms), I1.VCV1MORIG_EVALFUNC_CV2{h} );
                            I.VCV2MPERM_S(TsInd2,h,perms)       = cellmat_mergecols(I.VCV2MPERM_S(TsInd2,h,perms), num2cell(I1.DS_perm{h}(:,:,perms),2));
                        end
                    end
                    %I.VCV2MPERM_CV2(h,ll) = sum(I1.VCV1MPERM_CV2{h})/nperms(1);
                    clear Tx tmp V Ymodel modelTr modelTrL F Fkl dum 
                    %try, I = rmfield(I,'PCV1SUM'); end
                
                end
                I = nk_VisXHelper('accum', nM, nclass, decompfl, permfl, ix, jx, I, inp, ll, nperms(1), I1);  
                if any(permfl), 
                    for h=1:nclass
                        fprintf('\nCV2 [%g, %g]: Observed performance [ model #%g ]: %1.2f; P[CV2]: %1.3f', operm, ofold, h, I1.VCV1MORIG_EVALFUNC_CV2{h}, I.VCV2MPERM_CV2(h,ll) ); 
                    end
                end
                fprintf('\nSaving %s', oVISpath); save(oVISpath,'I1','sPs','operm','ofold');
                if saveparam, fprintf('\nSaving %s', OptModelPath); save(OptModelPath, 'MD', 'ofold','operm'); end
                if isfield(I1,'PCV1SUM'), PCV1SUMflag=true; else PCV1SUMflag = false; end
                clear Param I1 MD
                
            case 1

                vpth = deblank(vismat{f,d});

                if isempty(vpth) || ~exist(vpth,'file')
                    warning(['No valid VISdata-MAT detected for CV2 partition ' '[' num2str(f) ', ' num2str(d) ']!']);
                else
                    [~,vnam] = fileparts(vpth);
                    ind0(ll) = true; 
                    fprintf('\n\nLoading visualization data:');
                    fprintf('\n%s',vnam);
                    [I, I1] = nk_VisXHelper('accum', nM, nclass, decompfl, permfl, ix, jx, I, inp, ll, nperms(1), vpth);
                    if any(permfl)
                        for h=1:nclass
                             switch MODEFL
                                case 'classification'
                                    TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                                case 'regression'
                                    TsInd2 = CV.TestInd{f,d};
                            end
                            I.VCV2MORIG_S(TsInd2,h) = cellmat_mergecols(I.VCV2MORIG_S(TsInd2,h), num2cell(I1.DS{h},2));
                            for perms = 1:nperms(1), I.VCV2MPERM_S(TsInd2,h,perms) = cellmat_mergecols(I.VCV2MPERM_S(TsInd2,h,perms), num2cell(I1.DS_perm{h}(:,:,perms),2)); end
                        end
                    end
                    ol=ol+1;
                end
                if isfield(I1,'PCV1SUM'), PCV1SUMflag=true; else PCV1SUMflag = false; end

        end
        ll=ll+1; clear GDFEAT
    end
end


%%%%%%%%%%%%%%%%%%%%%%%% PERFORM IMAGING PROCEDURES %%%%%%%%%%%%%%%%%%%%%%%%%
if ~batchflag
    
    visdata = cell(1,nM);
    if any(permfl) 
        I.VCV2MODELP = nm_nanmedian(I.VCV2MPERM,2); 
        I.VCV2MODELP_STD = nm_nanstd(I.VCV2MPERM); 
        for h=1:nclass
            switch MODEFL
                case 'classification'
                    if numel(CV.class{1,1}{h}.groups) == 2
                        ind1 = inp.labels == CV.class{1,1}{h}.groups(1); f1 = ones(sum(ind1),1);
                        ind2 = inp.labels == CV.class{1,1}{h}.groups(2); f2 = -1*ones(sum(ind2),1);
                        labelh = zeros(numel(inp.labels),1);
                        labelh(ind1) = f1; labelh(ind2) = f2; %labelh(~labelh)=[];
                    else
                        labelh = zeros(size(inp.labels,1),1);
                        ind1 = inp.labels == CV.class{1,1}{h}.groups(1); 
                        labelh(ind1) = 1; labelh(~ind1,h) = -1;
                    end
                case 'regression'
                    labelh = inp.labels;
            end
            indempt = cellfun(@isempty,I.VCV2MORIG_S);
            Porig = cellfun(@nm_nanmedian,I.VCV2MORIG_S(~indempt(:,h),h)); 
            Lorig = labelh(~indempt(:,h));
            indnonnan = ~isnan(Porig); 
            if inp.targscale, IN.revertflag = true; IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; Porig = nk_PerfScaleObj(Porig, IN); end
            I.VCV2MORIG_EVALFUNC_GLOBAL(h) = feval(EVALFUNC, Lorig(indnonnan), Porig(indnonnan));
            fprintf('\nTesting observed model performance against permuted models using entire data: %g permutations\n', nperms(1))
            for perms = 1 : nperms(1)
                Pperm = cellfun(@nm_nanmedian, I.VCV2MPERM_S(~indempt(:,h),h,perms));
                if inp.targscale, Pperm = nk_PerfScaleObj(Pperm, IN); end
                I.VCV2MPERM_EVALFUNC_GLOBAL(h,perms) = feval(EVALFUNC, labelh(indnonnan), Pperm(indnonnan)); 
                crt = feval(compfun, I.VCV2MPERM_EVALFUNC_GLOBAL(h,perms), I.VCV2MORIG_EVALFUNC_GLOBAL(h));
                if ~crt, fprintf('*'); else fprintf('.'); end
                I.VCV2MPERM_GLOBAL(h,perms) = crt;
            end
        end
        if ~isempty(inp.extraL)
            I = nk_VisModels_ExtraLabels(I, inp, nperms, compfun);
        end
    end
    
    for n=1:nM
        
        % Number of classifiers / predictors
        [ D, datatype, brainmaski, badcoordsi, labeli, labelopi ] = getD(FUSION.flag, inp, n);
        
        if iscell(inp.VIS), nVIS = inp.VIS{n}; else nVIS = inp.VIS; end
        
        % Probability of feature selection across all CV2 * CV1 partitions
        for h=1:nclass
            NumPredDiv = repmat(I.VCV2NMODEL(h),size(I.VCV2SEL{h,n},1),1);
            if ~decompfl 
                if PCV1SUMflag && h==1, I.PCV2 = zeros(D,nclass); end
                I.PCV2(:,h) = I.PCV2SUM{h, n}' ./ NumPredDiv; 
                if size(I.VCV2PEARSON{h, n},2)>1
                    I.VCV2PEARSON_STD{h,n} = nm_nanstd(I.VCV2PEARSON{h, n},2);
                    I.VCV2PEARSON{h,n} = nm_nanmedian(I.VCV2PEARSON{h, n},2);
                    I.VCV2SPEARMAN_STD{h,n} = nm_nanstd(I.VCV2SPEARMAN{h, n},2);
                    I.VCV2SPEARMAN{h,n} = nm_nanmedian(I.VCV2SPEARMAN{h, n},2);
                    I.VCV2PEARSON_UNCORR_PVAL_STD{h,n} = nm_nanstd(I.VCV2PEARSON_UNCORR_PVAL{h, n},2);
                    I.VCV2PEARSON_UNCORR_PVAL{h,n} = nm_nanmedian(I.VCV2PEARSON_UNCORR_PVAL{h, n},2);
                    I.VCV2SPEARMAN_UNCORR_PVAL_STD{h,n} = nm_nanstd(I.VCV2SPEARMAN_UNCORR_PVAL{h, n},2);
                    I.VCV2SPEARMAN_UNCORR_PVAL{h,n} = nm_nanmedian(I.VCV2SPEARMAN_UNCORR_PVAL{h, n},2);
                    I.VCV2PEARSON_FDR_PVAL{h,n} = nm_nanmedian(I.VCV2PEARSON_FDR_PVAL{h, n},2);
                    I.VCV2SPEARMAN_FDR_PVAL{h,n} = nm_nanmedian(I.VCV2SPEARMAN_FDR_PVAL{h, n},2);
                    if linsvmfl
                        I.VCV2PVAL_ANALYTICAL{h,n} = nm_nanmedian(I.VCV2PVAL_ANALYTICAL{h, n},2);
                        I.VCV2PVAL_ANALYTICAL_FDR{h,n} = nm_nanmedian(I.VCV2PVAL_ANALYTICAL_FDR{h, n},2);
                    end
                end
            end

            if any(permfl),
                I.VCV2ZSCORE{h,n}   = nm_nansum(I.VCV2ZSCORE{h,n},2) ./ I.VCV2SEL{h,n};
                % Uncorrected p values
                Pvals = nm_nansum(I.VCV2PERM{h,n},2)./ol ;
                Pvals(Pvals==0) = realmin;
                %Pvals    = 1-normcdf(I.VCV2ZSCORE{h,n}); 
                I.VCV2PERM{h,n} = -log10(Pvals);
                % FDR-Corrected p values
                Pvals = nm_nansum(I.VCV2PERM_FDR{h,n},2)./ol ;
                Pvals(Pvals==0) = realmin;
                I.VCV2PERM_FDR_PVAL{h,n} = -log10(Pvals);
            end
            
            % Mean image across all CV2 * CV1 partitions
            I.VCV2{h,n} = nm_nansum(I.VCV2SUM{h,n},2) ./ I.VCV2SEL{h,n};
            
            % Compute Standard error 
            I.VCV2SUM2{h,n} = nm_nansum((I.VCV2SUM{h,n}.^2),2) ./ I.VCV2SEL{h,n};
            if size(I.VCV2SQ{h, n},2)>1
                I.VCV2SE{h,n}  = sqrt(abs(nm_nanmedian(I.VCV2SQ{h, n},2) - I.VCV2SUM2{h,n})./(I.VCV2SEL{h,n}-1));
                I.VCV2MEAN_CV1{h,n} = nm_nanmedian(I.VCV2MEAN{h,n},2);
                I.VCV2SE_CV1{h,n}   = nm_nanmedian(I.VCV2STD{h,n},2);
            else
                I.VCV2SE{h,n}  = sqrt(abs(I.VCV2SQ{h, n} - I.VCV2SUM2{h,n})./(I.VCV2SEL{h,n}-1));
                I.VCV2MEAN_CV1{h,n} = I.VCV2MEAN{h,n};
                I.VCV2SE_CV1{h,n}   = I.VCV2STD{h,n};
            end
            
            % Compute CV-ratio 
            I.VCV2rat{h,n} = I.VCV2{h,n} ./ I.VCV2SE{h,n};
            
             % Compute GrandMean metrics
            I.VCV2rat_CV1{h,n}  = I.VCV2MEAN_CV1{h,n}./I.VCV2SE_CV1{h,n};
            I.VCV2MEANthreshSE_CV1{h,n} = zeros(size(I.VCV2MEAN_CV1{h,n}),'single');
            indMEANgrSE = abs(I.VCV2MEAN_CV1{h,n}) > I.VCV2SE_CV1{h,n};
            I.VCV2MEANthreshSE_CV1{h,n}(indMEANgrSE) = I.VCV2MEAN_CV1{h,n}(indMEANgrSE);

            % Compute voxel selection probability using 95% confidence interval
            % method
            I.VCV2PROB{h,n} = (nm_nansum(I.VCV2PROB{h,n},2)/ol).*sign(I.VCV2MEAN_CV1{h,n});
        end

        if numel(inp.tF)>1 && datatype
            fprintf('\n\nWriting out images for Modality #%g',i)
        end
        
        % Now we have to differentiate between imaging and non-imaging
        % analyses. In the former case we write out data to the disk
        
        switch datatype
            % SPM-based NIFTI write-out
            case 1
                currdir = pwd;
                cd(inp.rootdir);
                for h=1:nclass % Loop through binary classifiedars (predictors)
                    % Generate filenames & save data
                    imgname = SAV.matname; 
                    suff = ['_NumPred-' num2str(I.VCV2NMODEL(h))];
                    varsuff = sprintf('_var%g',inp.tF(n));
                    switch MODEFL
                        case 'regression'
                                basename ='PredictVol';
                                suff = [multlabelstr suff varsuff '_ID' id];
                        case 'classification'
                                basename = 'DiscrimVol';
                                suff = [multlabelstr '_cl' num2str(h) suff varsuff '_ID' id];
                    end

                    % Save mean image
                    volnam = [basename '_Mean_' imgname suff ];
                    nk_WriteVol(I.VCV2{h,n},volnam,2, brainmaski,[], labeli, labelopi);

                    % Save CV ratio image
                    volnam = [basename '_CVratio_' imgname suff ];
                    nk_WriteVol(I.VCV2rat{h,n},volnam,2, brainmaski,[], labeli, labelopi);

                    % Save standard error image
                    volnam = [basename '_SE_' imgname suff ];
                    nk_WriteVol(I.VCV2SE{h,n},volnam,2, brainmaski,[], labeli, labelopi);

                    % Save grand mean image (CV1 level)
                    volnam = [basename '_Mean-GrM_' imgname suff ];
                    nk_WriteVol(I.VCV2MEAN_CV1{h,n},volnam,2, brainmaski,[], labeli, labelopi);

                    % Save grand mean standard error image (CV1 level)
                    volnam = [basename '_SE-GrM_' imgname suff ];
                    nk_WriteVol(I.VCV2SE_CV1{h,n},volnam,2, brainmaski,[], labeli, labelopi);

                    % Save grand mean thresholded by (mean > standard error) image (CV1 level)
                    volnam = [basename '_Mean-gr-SE-GrM_' imgname suff ];
                    nk_WriteVol(I.VCV2MEANthreshSE_CV1{h,n},volnam,2, brainmaski,[], labeli, labelopi);

                    % Save grand mean CV ratio image 
                    volnam = [basename '_CVratio-GrM_' imgname suff ];
                    nk_WriteVol(I.VCV2rat_CV1{h,n},volnam,2, brainmaski,[], labeli, labelopi);                   

                    % Save grand mean CV ratio image 
                    volnam = [basename '_Prob95CI-GrM_' imgname suff ];
                    nk_WriteVol(I.VCV2PROB{h,n},volnam,2, brainmaski,[], labeli, labelopi);  
                    
                    % Save grand mean Spearman image 
                    if exist('I.VCV2SPEARMAN','var')
                        volnam = [basename '_Spearman-GrM_' imgname suff ];
                        if ~isempty(I.VCV2SPEARMAN{h,n}), nk_WriteVol(I.VCV2SPEARMAN{h,n},volnam,2, brainmaski, [], labeli, labelopi);  end
                    end
                    % Save grand mean Spearman image 
                    if exist('I.VCV2PEARSON','var')
                        volnam = [basename '_Pearson-GrM_' imgname suff ];
                        if ~isempty(I.VCV2PEARSON{h,n}), nk_WriteVol(I.VCV2PEARSON{h,n},volnam,2, brainmaski,[], labeli, labelopi);  end
                    end
                    
                    if permfl 
                        
                        % Save permutation-based probability image
                        volnam = [basename '_PermProb_' imgname suff ];
                        nk_WriteVol(-log10(I.VCV2PERM{h,n}),volnam,2, brainmaski,[], labeli, labelopi);
                        volnam = [basename '_PermProbFDR_' imgname suff ];
                        nk_WriteVol(I.VCV2PERM_FDR_PVAL{h,n} ,volnam,2, brainmaski,[], labeli, labelopi);
                        
                        % Save permutation-based Zscore image 
                        volnam = [basename '_PermZ_' imgname suff ];
                        % Save permutation-based probability image
                        nk_WriteVol(I.VCV2ZSCORE{h,n},volnam,2, brainmaski,[], labeli, labelopi);
                    end

                    % Save grand mean image
                    if isfield(nVIS,'mean')
                        volnam = [basename '_GrandMean_thresh_' imgname suff ];
                        mI.GCV2_t = nk_Threshold(mI.GCV2, inp.VIS{n}.mean);
                        nk_WriteVol(mI.GCV2_t(:),volnam,2, brainmaski,[], labeli, labelopi);

                        % Save grand STD image
                        volnam = [basename '_GrandSTD_thresh_' imgname suff ];
                        sdI.GCV2_t = nk_Threshold(sdI.GCV2, inp.VIS{n}.se);
                        nk_WriteVol(sdI.GCV2_t(:),volnam,2, brainmaski,[], labeli, labelopi);

                        % Save grand Z score image
                        volnam = [basename '_GrandZ_thresh_' imgname suff ];
                        ind = ( mI.GCV2_t ~= 0 & sdI.GCV2_t > 0 );
                        normI.GCV2 = zeros(size(zI.GCV2));
                        normI.GCV2(ind) = zI.GCV2(ind);
                        nk_WriteVol(normI.GCV2(:),volnam,2, brainmaski,[], labeli, labelopi);

                        % Save grand probability image
                        indCV2MEANgrCV2SE = abs(I.GCV2SUMh) > repmat(seI.GCV2,1,ol);
                        GrandProb95CI = sum(indCV2MEANgrCV2SE,2)./ol;
                        volnam = [basename '_GrandProb95CI_' imgname suff ];
                        nk_WriteVol(GrandProb95CI(:),volnam,2, brainmaski,[], labeli, labelopi);
                    end

                    clear mI.GCV2_t sdI.GCV2_t ind normI.GCV2 zI.GCV2 indCV2MEANgrCV2SE GrandProb95CI

                    if isfield(nVIS,'thresh') && ~isempty(nVIS.thresh) && nVIS.thresh

                        I.VCV2_t = nk_Threshold(I.VCV2{h,n}, inp.VIS{n}.mean);
                        I.VCV2SE_t = nk_Threshold(I.VCV2SE{h,n}, inp.VIS{n}.se);
                        ind = ( I.VCV2_t ~= 0 & I.VCV2SE_t > 0 );

                        I.VCV2rat_t = zeros(size(I.VCV2_t));
                        I.VCV2rat_t(ind) = I.VCV2_t(ind) ./ I.VCV2SE_t(ind);

                        % Save thresholded mean image
                        threshstr = num2str(inp.VIS{n}.mean.val);
                        threshstr = regexprep(threshstr,' ','-');
                        volnam = [basename '_Mean_thresh-' threshstr '_' imgname suff ];
                        nk_WriteVol(I.VCV2_t,volnam,2, brainmaski,[], labeli, labelopi);

                        % Save thresholded CV ratio image
                        volnam = [basename '_CVratio_thresh_' imgname suff ];
                        nk_WriteVol(I.VCV2rat_t,volnam,2, brainmaski,[], labeli, labelopi);

                        % Save standard error image
                        threshstr = num2str(inp.VIS{n}.se.val);
                        threshstr = regexprep(threshstr,' ','-');
                        volnam = [basename '_SE_thresh-' threshstr '_' imgname suff ];
                        nk_WriteVol(I.VCV2SE_t,volnam,2, brainmaski,[], labeli, labelopi);
                    end 
                end
                cd(currdir);
            case 2 %Freesurfer write-out
        end 
        
         %% Build output structure
        visdata{n}.params.dimvecx      = [1 D];
        visdata{n}.params.varind       = inp.tF(n);
        visdata{n}.params.visflag      = datatype;
        visdata{n}.params.brainmask    = brainmaski;
        visdata{n}.params.badcoords    = badcoordsi;
        visdata{n}.params.I.numCV2part   = ll-1;
        visdata{n}.params.NumPred      = I.VCV2NMODEL(h);
        
        if ~isempty(featnames) && ~isempty(featnames{n});
            visdata{n}.params.features = featnames{n};
        else
            visdata{n}.params.features = cellstr(num2str((1:D)'));
        end
        visdata{n}.params.nfeats       = numel(visdata{n}.params.features);
        visdata{n}.MEAN                = I.VCV2(:,n);
        visdata{n}.SE                  = I.VCV2SE(:,n);
        visdata{n}.CVRatio             = I.VCV2rat(:,n);
        if isfield(I,'PCV2'), 
            visdata{n}.FeatProb        = {I.PCV2}; end
        visdata{n}.MEAN_CV2            = I.VCV2MEAN_CV1(:,n);
        visdata{n}.SE_CV2              = I.VCV2SE_CV1(:,n);
        visdata{n}.CVRatio_CV2         = I.VCV2rat_CV1(:,n);
        visdata{n}.Prob_CV2            = I.VCV2PROB(:,n);
        if ~decompfl(n)
            visdata{n}.Pearson_CV2              = I.VCV2PEARSON(:,n);
            visdata{n}.Spearman_CV2             = I.VCV2SPEARMAN(:,n);
            visdata{n}.Pearson_CV2_p_uncorr     = I.VCV2PEARSON_UNCORR_PVAL(:,n);
            visdata{n}.Spearman_CV2_p_uncorr    = I.VCV2SPEARMAN_UNCORR_PVAL(:,n);
            visdata{n}.Pearson_CV2_p_fdr        = I.VCV2PEARSON_FDR_PVAL(:,n);
            visdata{n}.Spearman_CV2_p_fdr       = I.VCV2SPEARMAN_FDR_PVAL(:,n);
            if isfield(I,'VCV2PEARSON_STD')
                visdata{n}.Pearson_CV2_STD          = I.VCV2PEARSON_STD(:,n);
                visdata{n}.Spearman_CV2_STD         = I.VCV2SPEARMAN_STD(:,n);
                visdata{n}.Pearson_CV2_p_uncorr_STD = I.VCV2PEARSON_UNCORR_PVAL_STD(:,n);
                visdata{n}.Spearman_CV2_p_uncorr_STD= I.VCV2SPEARMAN_UNCORR_PVAL_STD(:,n);
            end
            if linsvmfl
                visdata{n}.Analytical_p = I.VCV2PVAL_ANALYTICAL(:,n);
                visdata{n}.Analyitcal_p_fdr = I.VCV2PVAL_ANALYTICAL_FDR(:,n);
            end
        end
        if any(permfl)
            visdata{n}.PermProb_CV2             = I.VCV2PERM(:,n);
            visdata{n}.PermProb_CV2_FDR         = I.VCV2PERM_FDR(:,n); 
            visdata{n}.PermZ_CV2                = I.VCV2ZSCORE(:,n);
            visdata{n}.PermModel                = I.VCV2MODELP;
            visdata{n}.PermModel_std            = I.VCV2MODELP_STD;
            visdata{n}.PermModel_CV2            = I.VCV2MPERM_CV2;
            visdata{n}.PermModel_Eval_Global    = I.VCV2MPERM_GLOBAL;
            visdata{n}.PermModel_Crit_Global    = I.VCV2MPERM_EVALFUNC_GLOBAL;
            visdata{n}.ObsModel_Eval_Global     = I.VCV2MORIG_EVALFUNC_GLOBAL;
        end
        if ~isempty(inp.extraL)
            visdata{n}.ExtraL = I.EXTRA_L;
        end
    end
end
