%==========================================================================
%FORMAT [IN, OUT] = FoldPerm(IN, OUT, strout, RFEflag, FullPartflag, ...
%                            kxFlag, LoopParam)
%==========================================================================
%MAIN TRAINING MODULE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 07/2011

function [IN, OUT] = FoldPerm2(IN, OUT, strout, RFEflag, FullPartflag, RetrainImmediate, kxFlag, LoopParam)

global VERBOSE CV MODEFL RAND MULTILABEL

RF = []; 

%% Setup cell array containers
if ~exist('OUT', 'var') || isempty(OUT)
    
    OUT.tr      = cell( IN.nperms, IN.nfolds, IN.nclass ); % Training performance
    OUT.mdl     = cell( IN.nperms, IN.nfolds, IN.nclass ); % Model structures
    OUT.Trtargs = cell( IN.nperms, IN.nfolds, IN.nclass ); % Predicted target labels for CV1 training data
    OUT.Trdecs  = cell( IN.nperms, IN.nfolds, IN.nclass ); % Decision values / probabilities for CV1 training data
    OUT.kxVec   = cell( IN.nperms, IN.nfolds, IN.nclass ); % Subspace Stepping
    OUT.ts      = cell( IN.nperms, IN.nfolds, IN.nclass ); % Test performance
    OUT.rf      = cell( IN.nperms, IN.nfolds, IN.nclass ); % Recursive feature elimination structures (unused)
    OUT.w2      = cell( IN.nperms, IN.nfolds, IN.nclass ); % |w| for linear SVM kernels
    OUT.Md      = cell( IN.nperms, IN.nfolds, IN.nclass ); % distance to hyperplane for linear SVM kernels
    OUT.Mm      = cell( IN.nperms, IN.nfolds, IN.nclass ); % normalized margin for linear SVM kernels
    OUT.CVtargs = cell( IN.nperms, IN.nfolds, IN.nclass ); % Predicted target labels for CV1 test data
    OUT.CVdecs  = cell( IN.nperms, IN.nfolds, IN.nclass ); % Decision values / probabilities for CV1 test data
    OUT.featnum = cell( IN.nperms, IN.nfolds, IN.nclass, IN.nvar ); % Number of feature
    
end
nc = size(IN.F,3); % number of binary comparisons in IN.F (features)

if exist('LoopParam','var') && ~isempty(LoopParam)
    PermVec = LoopParam.PermVec;
    FoldVec = LoopParam.FoldVec;
    ClassVec = LoopParam.ClassVec;
else
    PermVec = 1:IN.nperms;
    FoldVec = 1:IN.nfolds; 
    ClassVec = 1:IN.nclass;
end

PermNum = numel(PermVec);
FoldNum = numel(FoldVec);
ClassNum = numel(ClassVec);
kx = 1;

for ii=1:PermNum % Loop through CV1 permutations

    for jj=1:FoldNum % Loop through CV1 folds
        
        i = PermVec(ii); j= FoldVec(jj);
        
        for ccurclass=1:ClassNum % Loop through dichotomizers
            
            curclass = ClassVec(ccurclass);
            
            tFea = zeros(1, IN.nvar); 
            Fk = cell(1, IN.nvar); 
            tCV = Fk; 
            tTr = Fk; 
            modelTr = Fk;
            CVInd = IN.Y.CVInd{i,j}{curclass};
            TrInd = IN.Y.TrInd{i,j}{curclass};
            
            %% Loop through variates and extract multi-modal training data cell array
            for v=1:IN.nvar 
                
                % Get feature subspace size
                if nc > 1 % more than one binary classifier
                    try 
                        tFea(v) = size(IN.F{i,j,curclass,v},2); 
                    catch
                        tFea(v) =  size(IN.F{i,j,curclass,v},1); 
                    end
                    Fk{v} = IN.F{i,j,curclass,v}; 
                else
                    try % a binary classifier (a multi-group classifier, not impl.)
                        tFea(v) = size(IN.F{i,j,v},2); 
                    catch
                        tFea(v) =  size(IN.F{i,j,v},1); 
                    end
                    Fk{v} = IN.F{i,j,v}; 
                end
                
                % Determine training data & validation data
                if ~iscell(IN.Y.Tr{i,j,1})
                    modelTr{v} = IN.Y.Tr{i,j,v}(TrInd,:);
                    tTr{v} = IN.Y.Tr{i,j,v};
                    if FullPartflag, modelTr{v} = [modelTr{v}; IN.Y.CV{i,j,v}(CVInd,:) ]; end
                    tCV{v} = IN.Y.CV{i,j,v};
                else
                    modelTr{v} = IN.Y.Tr{i,j,v}{curclass}(TrInd,:);
                    tTr{v} = IN.Y.Tr{i,j,v}{curclass};
                    if FullPartflag, modelTr{v} = [modelTr{v}; IN.Y.CV{i,j,v}{curclass}(CVInd,:) ]; end
                    tCV{v} = IN.Y.CV{i,j,v}{curclass};
                end
            end
            
            %% Get samples sizes
            kSubjTr = size(tTr{1},1); 
            kSubjCV = size(tCV{1},1);
            
            %% Determine labels for learning process
            modelTrL = IN.Y.TrL{i,j}{curclass}(:,MULTILABEL.curdim);
            tTrL = zeros(size(tTr{1},1),1); 
            tTrL(TrInd) = IN.Y.TrL{i,j}{curclass}(:,MULTILABEL.curdim);
            tCVL = zeros(size(tCV{1},1),1); 
            tCVL(CVInd) = IN.Y.CVL{i,j}{curclass}(:,MULTILABEL.curdim);
            if FullPartflag, 
                modelTrL = [modelTrL; IN.Y.CVL{i,j}{curclass}(:,MULTILABEL.curdim)]; 
            end
            
            %% Define feature subspace counter ...
            % according to maximum
            % feature subspace size across variates. Compute subspace
            % stepping for loop.
            % Stepping for feature subspace search
            lFea = max(tFea); 
            
            if kxFlag, kx = ceil((lFea / 100) * kxFlag); end
            
            OUT.kxVec{i,j,curclass} = kx:kx:lFea;
            if isempty(OUT.kxVec{i,j,curclass})
                OUT.kxVec{i,j,curclass} = 1;
            else
                if OUT.kxVec{i,j,curclass}(end) < lFea, 
                    OUT.kxVec{i,j,curclass} = [ OUT.kxVec{i,j,curclass} lFea ]; 
                end; 
            end
            kFea = length(OUT.kxVec{i,j,curclass});
            cvts_fl = kFea ~= size(OUT.ts{i,j,curclass},1);
            %% Initialize arrays for current multi-dimensional cell pointer
            OUT.tr{i,j,curclass}         = zeros( kFea, 1 );
            if isempty(OUT.mdl{i,j,curclass})
                OUT.mdl{i,j,curclass}        = cell( kFea, 1 );
            end
            OUT.Trtargs{i,j,curclass}    = zeros( kSubjTr, kFea );
            OUT.Trdecs{i,j,curclass}     = zeros( kSubjTr, kFea );
            
            if ~FullPartflag || RetrainImmediate || cvts_fl
                OUT.ts{i,j,curclass}     = zeros( kFea,1 );
                OUT.CVtargs{i,j,curclass}= zeros( kSubjCV, kFea );
                OUT.CVdecs{i,j,curclass} = zeros( kSubjCV, kFea );
                OUT.rf{i,j,curclass}     = cell( kFea, 1 );
                OUT.w2{i,j,curclass}     = zeros( kFea, 1 );
                OUT.Md{i,j,curclass}     = zeros( kFea, 1 );
                OUT.Mm{i,j,curclass}     = zeros( kFea, 1 );
            end
            
            % Loop through every (or every kx, to save run time) feature subspace
            indkX = 1;
            
            for kT=1:kFea
                
                k = OUT.kxVec{i,j,curclass}(kT);
                nfeats = sum(any(Fk{1}(:,k),2));
                if VERBOSE, 
                    switch MODEFL
                        case 'classification'
                           if RAND.Decompose ~= 9
                                fprintf('\n%s => CV1 [%g, %g, %s, Suspace %g/%g (%g feats)]:', strout, ...
                                    i, j, CV.class{1,1}{curclass}.groupdesc, k, lFea, nfeats)
                           else
                               fprintf('\n%s => CV1 [%g, %g, Multi-Group, Subspace %g/%g (%g feats)]:', strout, ...
                                    i, j, k, lFea, nfeats)
                           end
                        case 'regression'
                           fprintf('\n%s => CV1 [%g, %g, Regression, Subspace %g/%g (%g feats)]:', strout, ...
                                    i, j, k, lFea, nfeats)
                    end
                end
                
                %% Extract training and test data according to current feature subspace mask Fk
                [Ymodel, Fx] = nk_ExtractFeatures(modelTr, Fk, [], k);
                Ytrain = nk_ExtractFeatures(tTr, Fk, [], k);
                Ytest  = nk_ExtractFeatures(tCV, Fk, [], k);
                
                %% Train model(s)
                switch RFEflag
                    
                   case 0 % no RFE
                       
                        % Train algorithm 
                        [~, model]= nk_GetParam2(Ymodel, modelTrL, IN.Ps{curclass}, 1);
                        
                    case 1 % some RFE type
                       
                        if IN.nvar < 2 % Univariate case
                            OUT.featout{i,j,curclass} = zeros(size(IN.F{i,j,curclass},1),kFea); 
                        else % Multiple variate case
                            OUT.featout = cell(IN.nvar,1);
                            for v=1:IN.nvar 
                                OUT.featout{v}{i,j,curclass} = zeros(size(IN.F{i,j,curclass,v},1),tFea(v)); 
                            end
                        end
                        
                        % NOTE: nk_RFE currently supports only univariate data
                        [RF, model]  = nk_MLOptimizer_Wrapper(Ymodel, modelTrL, Ytest, tCVL, IN.Ps{curclass}, []);
                        if isempty(RF.FeatureIndex), error('Wrapper did not returned any features! Relax your wrapper settings'); end
                        if RF.found
                            % Transfer new subspace features to SubSets 
                            indFx                       = find(Fx); 
                            indFx                       = indFx(RF.FeatureIndex);
                            indFx0                      = zeros(size(Fx));  
                            indFx0(indFx)               = IN.F{i,j,curclass,v}(indFx,k);
                            IN.F{i,j,curclass,v}(:,k)   = indFx0;
                            % Extract RFE subspace from Ytrain and Ytest
                            Ytrain                      = Ytrain(:,RF.FeatureIndex);
                            Ytest                       = Ytest(:,RF.FeatureIndex);
                            Ymodel                      = Ymodel(:,RF.FeatureIndex);
                        end
                end
                
                %% Apply trained algorithm to CV1 test data
                [OUT.tr{i,j,curclass}(indkX), ...
                    OUT.Trtargs{i,j,curclass}(:,indkX) , ...
                    OUT.Trdecs{i,j,curclass}(:,indkX) ] = nk_GetTestPerf(Ytrain, tTrL, [], model, Ymodel);
                if isnan(OUT.tr{i,j,curclass}(indkX))
                    warning('Non-finite performance measures found in CV1 training data')
                end
                if VERBOSE, fprintf('\tTr = %1.2f', OUT.tr{i,j,curclass}(indkX)); end
                if ~FullPartflag || RetrainImmediate || cvts_fl
                    [OUT.ts{i,j,curclass}(indkX), ...
                        OUT.CVtargs{i,j,curclass}(:,indkX), ...
                        OUT.CVdecs{i,j,curclass}(:,indkX)] = nk_GetTestPerf(Ytest, tCVL, [], model, Ymodel);
                    if isnan(OUT.ts{i,j,curclass}(indkX))
                        warning('Non-finite performance measures found in CV1 test data')
                    end
                    if VERBOSE, fprintf(', CV = %1.2f',OUT.ts{i,j,curclass}(indkX)); end
                end
                OUT.mdl{i,j,curclass}{indkX} = model;
                if ~isempty(RF) && RF.found
                    OUT.featout{i,j,curclass}(:,indkX) = indFx0;
                    OUT.featnum{i,j,curclass}(indkX) = sum(OUT.featout{i,j,curclass}(:,indkX)~=0);
                else
                    for v=1:IN.nvar 
                        OUT.featnum{i,j,curclass,v}(indkX) = sum(Fk{v}(:,indkX)~=0);
                    end
                end
                indkX = indkX+1;
            end

            %% Convert to single to save disk space
            OUT.Trtargs{i,j,curclass}    = single(OUT.Trtargs{i,j,curclass}); 
            OUT.Trdecs{i,j, curclass}    = single(OUT.Trdecs{i,j,curclass});
            OUT.tr{i,j, curclass}        = single(OUT.tr{i,j,curclass}); 
            if ~FullPartflag || RetrainImmediate || cvts_fl
                OUT.CVtargs{i,j, curclass}   = single(OUT.CVtargs{i,j,curclass}); 
                OUT.CVdecs{i,j, curclass}    = single(OUT.CVdecs{i,j,curclass});
                OUT.ts{i,j, curclass}        = single(OUT.ts{i,j,curclass}); 
            end
        end
    end
end

%if detrendfl
%
%end

end