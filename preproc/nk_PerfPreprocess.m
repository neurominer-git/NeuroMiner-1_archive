function [tY, Pnt, paramfl, tYocv] = nk_PerfPreprocess(Y, inp, labels, paramfl, Yocv, Cocv)
% =========================================================================
% function [tY, Param] = nk_PerfPreprocess(Y, i, j, kbin, labels, ...
%                                                      paramfl, Yocv, Cocv)
% =========================================================================
%
% INPUTS:
% -------
% Y         = data matrix as [m x n] double, m = samples, n = dimensions
% inp 		= înput parameter structure
% labels    = n x 1 label vector 
% paramfl   = script-controlling parameters 
% Yocv      = Independent test data
% Cocv      = Calibration data
% 
% OUTPUTS:
% --------
% tY        = the preprocessed data
% Param     = computed preprocessing parameters
% paramfl   = modified running parameters
% tYocv     = the preprocessed independent test data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2017

global PREPROC MODEFL MULTI CV RAND VERBOSE TEMPL SAV

% Initialize runtime variables
i       = inp.f; 
j       = inp.d; 
kbin    = inp.nclass;
[iy,jy] = size(CV(1).cvin{i,j}.TrainInd); 
TrInd   = CV(1).TrainInd{i,j}; 
TsInd   = CV(1).TestInd{i,j}; 
TsI     = cell(kbin,1);
if ~exist('paramfl','var'), paramfl.use_exist = false; end;
cv2flag = false; if isfield(PREPROC,'CV2flag') && (PREPROC.CV2flag - 1) == 1; cv2flag = true; end
tYocv   = []; if exist('Yocv','var') && ~isempty(Yocv), tYocv.Ts = cell(iy,jy); else, Yocv = []; end
if ~exist('Cocv','var'), Cocv = []; end

% Set binary/regression or multi-group processing mode
if iscell(PREPROC)
    BINMOD = PREPROC{1}.BINMOD;
else
    BINMOD = PREPROC.BINMOD;
end
if BINMOD || strcmp(MODEFL,'regression')
    ukbin = kbin;   if VERBOSE; fprintf('\nProcessing Mode: binary / regression preprocessing'); end
else
    ukbin = 1;      
    if VERBOSE; fprintf('\nProcessing Mode: multi-group preprocessing'); end
end

if VERBOSE
    if iscell(Y),
        fprintf('\nMultiple shelfs of input data detected')
        for ii=1:size(Y,1)
            fprintf('\nShelf [ %2g ]: Original dimensionality: %g', ii, size(Y{ii},2)); 
        end
    else
        fprintf('\nOriginal dimensionality: %g', size(Y,2)); 
    end    
end

% Set up data containers
tY.Tr = cell(iy,jy);
tY.CV = cell(iy,jy);
tY.Ts = cell(iy,jy);
    
% Labels & Indices
tY.TrL = cell(iy,jy);
tY.CVL = cell(iy,jy);
tY.TrInd = cell(iy,jy);
tY.CVInd = cell(iy,jy);

% Eventually, apply spatial operations to image data
% (This function will be extended beyond smoothing ops on nifti data)
[sY, sYocv, sCocv, inp, optfl, ukbin, uBINMOD, BINMOD] = ...
    nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, paramfl, BINMOD, kbin, ukbin);

if ~BINMOD && isfield(paramfl,'PXopt') && numel(paramfl.PXopt)>1,
    uBINMOD = 0; ukbin = kbin;
    if VERBOSE; fprintf('\nProcessing Mode: multi-group preprocessing, but no multi-group classifier requested'); end
end

% Generate template parameters (e.g. for Procrustes rotation)
if isfield(paramfl,'templateflag') && paramfl.templateflag 
    TEMPL = nk_GenTemplParam(PREPROC, CV, MODEFL, RAND, sY, inp, kbin);
else
    clear TEMPL
end

% The order of statements here is critical: TsI should receive the original
% TsInd before it is modified below!
if ~isempty(MULTI) && MULTI.flag, tY.mTsL = labels(TsInd,:); end
for u=1:ukbin, TsI{u} = TsInd; end

for u=1:kbin
    switch MODEFL
        case 'regression'
            TsL     = labels(TsInd,:);
            TsInd   = true(size(TsL,1),1);
        case 'classification'
            if RAND.Decompose == 9
                TsL     = labels(TsInd,:);
                TsInd   = true(size(TsL,1),1);
            else
                TsL     = CV.classnew{i,j}{u}.label;
                TsInd   = CV.classnew{i,j}{u}.ind;
            end
    end
    tY.TsL{u} = TsL;
    tY.TsInd{u} = TsInd;
end

if ~exist('paramfl','var') || ~paramfl.found || ~paramfl.use_exist
    Pnts = struct('data_ind', [], ...
                  'train_ind', [], ...
                  'nP', [], ...
                  'nA', [], ...
                  'TrainedParam', []);
    Pnt = repmat(Pnts,iy,jy,ukbin);
else
    Pnt = paramfl.Param;
    paramfl.Param = [];
end

for k=1:iy % Inner permutation loop

    for l=1:jy % Inner CV fold loop
        
        for u=1:ukbin % Binary comparison loop depending on PREPROC.BINMOD & FBINMOD
            
            %Check whether optimised parameter space exists
            PREPROC = nk_SetParamChain(paramfl, u, PREPROC);

            if optfl
                usY = sY{u}; 
                if ~isempty(Yocv), usYocv = sYocv{u}; end
                if ~isempty(Cocv), usCocv = sCocv{u}; end
            else
                usY = sY;
                if ~isempty(Yocv), usYocv = sYocv; end
                if ~isempty(Cocv), usCocv = sCocv; end
            end
            paramfl.P{u} = nk_ReturnParamChain(PREPROC, 1); 
            
            %% Push data into partition / folds
            % Define pointers to data
            if ~uBINMOD || strcmp(MODEFL,'regression')
                % Multi-group label: 1, 2 ,3, ...
                TrX = TrInd(CV.cvin{i,j}.TrainInd{k,l});
            else
                % Binary label
                TrX = TrInd(CV.class{i,j}{u}.TrainInd{k,l});
            end
            TrI = TrInd(CV.cvin{i,j}.TrainInd{k,l});
            CVI = TrInd(CV.cvin{i,j}.TestInd{k,l});
            TrL = labels(TrI,:);
            CVL = labels(CVI,:);
            
            if size(TrI,2)>2
                TrI = TrI';
                CVI = CVI';
                TsI{u} = TsI{u}';
            end

            % Assign data 
            if iscell(usY)
                vTr = cell(numel(usY),1); if ~isempty(Yocv), vTs = cell(numel(usY),4); else vTs = cell(numel(usY),3); end
                for pu = 1:numel(usY);
                    % Training & CV data
                    vTr{pu} = usY{pu}(TrX,:); vTs{pu,1} = usY{pu}(TrI,:); vTs{pu,2} = usY{pu}(CVI,:); vTs{pu,3} = usY{pu}(TsI{u},:); 
                    % Independent test data
                    if ~isempty(Yocv), vTs{pu,4} = usYocv{pu}; end
                    % Calibration data
                    if ~isempty(Cocv), InputParam.C{pu} = usCocv{pu}; end
                           
                    % Remove cases which are completely NaN
                    [vTr{pu},~, SrcParam.iTrX] = nk_ManageNanCases(vTr{pu}, TrL);
                    [vTs{pu,1}, TrL, SrcParam.iTr] = nk_ManageNanCases(vTs{pu,1}, TrL);
                    [vTs{pu,2}, CVL, SrcParam.iCV] = nk_ManageNanCases(vTs{pu,2}, CVL);
                    [vTs{pu,3}, ~, SrcParam.iTs] = nk_ManageNanCases(vTs{pu,3});
                end
            else
                % Training & CV data
                vTr = usY(TrX,:); vTs{1} = usY(TrI,:); vTs{2} = usY(CVI,:); vTs{3} = usY(TsI{u},:);
                
                % Remove cases which are completely NaN
                [vTr,~,SrcParam.iTrX] = nk_ManageNanCases(vTr, TrL);
                [vTs{1}, TrL, SrcParam.iTr] = nk_ManageNanCases(vTs{1}, TrL);
                [vTs{2}, CVL, SrcParam.iCV] = nk_ManageNanCases(vTs{2}, CVL);
                [vTs{3}, ~, SrcParam.iTs] = nk_ManageNanCases(vTs{3});
                
                % Independent test data
                if ~isempty(Yocv), vTs{4} = usYocv; end
                % Calibration data
                if ~isempty(Cocv), InputParam.C = usCocv; end
            end
            if isfield(inp,'Yw'), InputParam.Yw = inp.Yw; end 
            InputParam.Tr = vTr; InputParam.Ts = vTs;

            if VERBOSE; fprintf('\n-----------------------------------------------------------------------------------'); 
                switch BINMOD
                    case {1,3,4}
                        if strcmp(MODEFL,'regression') || RAND.Decompose == 9
                            fprintf('\nWorking on data partition: CV2 [%g, %g], CV1 [%g, %g]', ...
                                i, j, k, l)
                        else
                            fprintf('\nWorking on data partition: CV2 [%g, %g], CV1 [%g, %g, %s]', ...
                                i, j, k, l, CV.class{i,j}{u}.groupdesc)
                        end
                    case {0,2}
                            fprintf('\nWorking on data partition: CV2 [%g, %g], CV1 [%g, %g]', ...
                                i, j, k, l)
                end
            end
            
            %% Generate SrcParam structure
            SrcParam.TrX                = TrX;
            SrcParam.TrI                = TrI;
            SrcParam.CVI                = CVI;
            SrcParam.TsI                = TsI{u};
            SrcParam.u                  = u;
            SrcParam.binmult            = BINMOD;
            SrcParam.CV1perm            = k;
            SrcParam.CV1fold            = l;
            SrcParam.covars             = inp.covars;
            SrcParam.covars_oocv        = inp.covars_oocv;
            
            %SrcParam.modalparams        = paramfl.PV;
            %SrcParam.modalparams        = rmfield(SrcParam.paramfl,'Y');
            
            switch MODEFL
                case 'classification'
                    if RAND.Decompose ~=9
                        SrcParam.BinaryTrainLabel   = CV.class{i,j}{u}.TrainLabel{k,l};
                        SrcParam.BinaryCVLabel      = CV.class{i,j}{u}.TestLabel{k,l};
                    end
                    SrcParam.MultiTrainLabel    = TrL;
                    SrcParam.MultiCVLabel       = CVL;
                case 'regression'
                    SrcParam.TrainLabel         = TrL;
                    SrcParam.CVLabel            = CVL;
            end

            %% Generate and execute for given CV1 partition preprocessing sequence

            [InputParam, oTrainedParam, SrcParam] = nk_GenPreprocSequence(InputParam, PREPROC, SrcParam, Pnt(k, l, u).TrainedParam);
            
            % Check whether we have imputed labels
            if isfield(SrcParam,'TrL_imputed'), TrL = SrcParam.TrL_imputed; end
            
            %% Write preprocessed data to mapY structure
            if isfield(paramfl,'PREPROC') && isfield(paramfl,'PXfull') && ~isempty(paramfl.PXopt{u})
                % Here an optimized parameter space exists that has
                % been used to limit the computation to the unique
                % parameter cominations. We create and store the pointers
                % and used them later on to retrieve the right preproc
                % version of the data and the respective preproc params
                [ Pnt(k,l,u).data_ind, ...
                  Pnt(k,l,u).train_ind, ...
                  Pnt(k,l,u).nP, ...
                  Pnt(k,l,u).nA] = nk_ParamReplicator(paramfl.P{u}, paramfl.PXopt{u}, paramfl.PREPROC, oTrainedParam);
            end
            trd = InputParam.Ts(:,1); 
            cvd = InputParam.Ts(:,2); 
            tsd = InputParam.Ts(:,3);
            if ~isempty(Yocv), ocv = InputParam.Ts(:,4); end
            if  any(cellfun(@isempty,trd)) ||  any(cellfun(@isempty,cvd)) ||  any(cellfun(@isempty,tsd)),
                error('Empty training/test/validation data matrices return by preprocessing pipeline. Check your settings')
            end
            clear InputParam
            
            TrL = labels(TrI,:); CVL = labels(CVI,:);
            switch BINMOD

                case 0 % Multi-group mode both in FBINMOD and PREPROC.BINMOD
                    if ukbin > 1
                        tY.Tr{k,l}{u} = nk_ManageNanCases(trd, [], SrcParam.iTr);
                        tY.CV{k,l}{u} = nk_ManageNanCases(cvd, [], SrcParam.iCV);
                        tY.Ts{k,l}{u} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                        if ~isempty(Yocv), tYocv.Ts{k,l}{u} = ocv; end
                    else
                        tY.Tr{k,l} = nk_ManageNanCases(trd, [], SrcParam.iTr);
                        tY.CV{k,l} = nk_ManageNanCases(cvd, [], SrcParam.iCV);
                        tY.Ts{k,l} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                        if ~isempty(Yocv), tYocv.Ts{k,l} = ocv; end
                    end
                    % Write dichotomization labels to CV1 partition
                    for zu=1:kbin % Binary loop depending on the # of binary comparisons
                        % Generate logical indices
                        if isfield(CV,'class') && length(CV.class{i,j}{zu}.groups) == 2
                            % One-vs-One
                            indtr = ( TrL == CV.class{i,j}{zu}.groups(1) | TrL == CV.class{i,j}{zu}.groups(2) ) | isnan(TrL);   
                            indcv = ( CVL == CV.class{i,j}{zu}.groups(1) | CVL == CV.class{i,j}{zu}.groups(2) ) | isnan(CVL);
                        else
                            % One-vs-All
                            indtr = TrL~=0;   
                            indcv = CVL~=0;
                        end
                        % Write indices
                        tY.TrInd{k,l}{zu} = indtr;
                        tY.CVInd{k,l}{zu} = indcv;
                        % Write labels to CV1 partition
                        if isfield(CV,'class')
                            tY.TrL{k,l}{zu} = CV.class{i,j}{zu}.TrainLabel{k,l};
                            tY.CVL{k,l}{zu} = CV.class{i,j}{zu}.TestLabel{k,l};	
                        else
                            tY.TrL{k,l}{zu} = labels(CV.TrainInd{i,j}(CV.cvin{i,j}.TrainInd{k,l}));
                            tY.CVL{k,l}{zu} = labels(CV.TrainInd{i,j}(CV.cvin{i,j}.TestInd{k,l}));
                        end
                    end

                case 1   

                    % Write data to CV1 partition
                    tY.Tr{k,l}{u} = nk_ManageNanCases(trd, [], SrcParam.iTr); 
                    tY.CV{k,l}{u} = nk_ManageNanCases(cvd, [], SrcParam.iCV); 
                    tY.Ts{k,l}{u} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                    if ~isempty(Yocv), tYocv.Ts{k,l}{u} = ocv; end
                    
                    if ~strcmp(MODEFL,'regression') && length(CV.class{i,j}{u}.groups) == 2
                        indtr = ( TrL == CV.class{i,j}{u}.groups(1) | TrL == CV.class{i,j}{u}.groups(2) ) | isnan(TrL);    
                        indcv = ( CVL == CV.class{i,j}{u}.groups(1) | CVL == CV.class{i,j}{u}.groups(2) ) | isnan(CVL);
                    else
                        indtr = true(size(TrI,1),1);   
                        indcv = true(size(CVI,1),1);
                    end
                    % Write indices
                    tY.TrInd{k,l}{u} = indtr;
                    tY.CVInd{k,l}{u} = indcv;

                    switch MODEFL
                        case 'regression' 
                            tY.TrL{k,l}{u} = TrL;
                            tY.CVL{k,l}{u} = CVL;
                        case 'classification'
                            if RAND.Decompose ~=9
                                % Write labels to CV1 partition
                                tY.TrL{k,l}{u} = CV.class{i,j}{u}.TrainLabel{k,l};	
                                tY.CVL{k,l}{u} = CV.class{i,j}{u}.TestLabel{k,l};	
                            else
                                tY.TrL{k,l}{u} = TrL;
                                tY.CVL{k,l}{u} = labels(CV.TrainInd{i,j}(CV.cvin{i,j}.TestInd{k,l}),:);
                            end
                    end

            end
            if ~isempty(MULTI) && MULTI.flag, tY.mTrL{k,l} = TrL; tY.mCVL{k,l} = CVL; end
            
            % save parameters
            if paramfl.write || cv2flag, Pnt(k,l,u).TrainedParam = oTrainedParam; end
            clear TrainedParam SrcParam
        end
    end
end

