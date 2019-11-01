function [ mW, mP, mR, mSR, W, mPA ]= nk_VisXWeight2(inp, MD, Y, L, varind, P, F, VI, decompfl, procfl, Fadd)
% ================================================================================
% [mW, mP, W] = nk_VisXWeight2(MD, Y, L, varind, P, F, VI, decompfl, procfl, Fadd)
% ================================================================================
% Core visualization module that retrieves a weight vector and maps it back
% to the input space of features by reversing processing steps as much as
% possible and meaningful
%
% Inputs:
% ------
% inp :
% MD :          Model
% Y :           Features training data
% L :           Labels training data
% varind :      Selected modality
% P :           Preprocessing parameters for given training data partition
% F :           Mask of selected features across modality
% VI :        
% decompfl :    Presence of dimensionality reduction method
% procfl :
% Fadd :
%
% Outputs:
% --------
% mW :          (Back-projected) weight map
% mP :
% mR :
% mSR :
% W :           weight vector in processed feature space
% mPA :
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2019

global SVM %TEMPL

if ~exist('procfl','var') || isempty(procfl), procfl = true; end
if ~exist('Fadd','var') || isempty(Fadd), Fadd = true(size(F)); end
mPA = []; PA = []; AnalP=[]; warning off
% Check what type of kernel has been used
switch SVM.prog
    case {'LIBSVM','LIBLIN'}
        switch SVM.kernel.kernstr
            case{' -t 0',' -t 4',' -t 5','lin', 'linear', 'none'} 
                %%%%%%%%%%%%% USE WEIGHTS OF MODEL %%%%%%%%%%%%
                xV = nk_GetPrimalW(MD); % Get weight vector over feature space
                if ~decompfl(1), 
                    % Remove cases which are completely NaN
                    [Yn, Ln] = nk_ManageNanCases(Y, L); 
                    AnalP = compute_analytical_pvals(Ln,Yn,xV'); 
                end
            otherwise % non-linear
                %%%%%%%%%% COMPUTE MIN. DIFF. VECTORS %%%%%%%%% 
                xV = nk_VisSV(MD, Y, L);
        end
    case 'MEXELM'
        xV = MD.outW;
    case 'GLMFIT'
        xV = MD.beta(2:end);
    case 'matLRN'
        % Check whether addBias == 1 and whether algo was kernalized
        if isfield(MD,'addBias') && MD.addbias, xV = MD.w(2:end); else, xV = MD.w; end
        if isfield(MD,'kernel') && MD.kernel, error('Unfortunately, the pre-image of the kernel weight vector cannot be computed with the current version NM.'); end
    case 'GLMNET'
        xV = MD.beta(:,end);
    case 'RNDFOR'
        xV = MD.importance;
    case 'MikRVM'
        % Use the Quality factor (=importance of feature for reducing the prediction error)
        % in the RVM model as weight vector in the subsequent computations
        xV = MD.D.Q_Factor';
    case 'DECTRE'
        xV = MD.predictorImportance;
    otherwise
        error(['Vizualisation for ' SVM.prog ' not supported yet. Please consult the NM Manual for more information']);
end
W       = zeros(size(F,1),1);
W(F)    = xV;
W(~Fadd)= 0;
if ~any(W), fprintf('?'); end

% Process analytical P values
if ~isempty(AnalP)
    PA = zeros(size(F,1),1);
    PA(F) = AnalP;
    PA(~Fadd) = 0;
end

Fu      = F & Fadd;
nM = numel(inp.PREPROC);
mM = numel(inp.tF);

if inp.norm == 1 && any(strcmp({'LIBSVM','LIBLIN'},SVM.prog))
    % The criterion for pattern significance in SVM is the weight vector normalized to
    % the margin width as described in Gaonkar et al. 2015, Medical Image Analysis
    % This however may result in overconservative P value estimation
    % because it flattens the null distribution
    W = W/(norm(W,2));
end

mW = cell(1,mM);
if ~isempty(PA), mPA = cell(1,mM); end

if procfl
    mP = cell(1,mM); mR = cell(1,mM); mSR = cell(1,mM);
else
    mP = []; mR = []; mSR = [];
end

% Loop through modalities
for n=1:nM
    
    if iscell(inp.PREPROC), 
        nPREPROC = inp.PREPROC{n};
        nPX      = P{n};
    else
        nPREPROC = inp.PREPROC;
        nPX      = P;
    end
    
    if nM > 1
        % Input spaces were concatenated before training the prediction algorithm 
        % (e.g. PCA was applied independently to each data space and the reduced feature spaces were then concatenated)
        % in NM we call this an intermediate fusion strategy and the
        % processing of the modality-specific subspaces data has to be done
        % independently (e.g. a separate PCA for each subspace)
        lVI = VI == n;  
    else
        % if early fusion a single preprocessing pipeline has been applied to
        % the concatenated feature spaces
        lVI = true(size(Fu)); 
    end
    
    lFuVI = Fu & lVI; fVI = find(lVI); nmP = Fu(fVI);
    nmW = zeros(numel(fVI),1); nmW(nmP) = W(lFuVI);
    if ~isempty(PA)
         nmPA = nan(numel(fVI),1); nmPA(nmP) = PA(lFuVI);
    end
    
    if procfl 
        
        if isfield(nPREPROC,'ACTPARAM'), nA = numel(nPREPROC.ACTPARAM); else, nA=0; end
        
        %%%%%%%%%%% PROJECT DATA BACK TO INPUT SPACE %%%%%%%%%%

        % Revert dimensionality reduction if previously used
        % Find Dimensionality reduction parameters
        reducedimfl = false; flg = false; 
        
        for a = nA:-1:1
            
            % Adjust pnt according to parameter combination in current
            % preprocessing step
            
            switch nPREPROC.ACTPARAM{a}.cmd

                case 'scale'
                    % Very frequently users opt to scale features after
                    % factorization/dimensionality reduction. In these cases
                    % scaling has to be reverted prior to back-projection
                    if ~reducedimfl && decompfl(n),
                        IN = nPREPROC.ACTPARAM{a}.SCALE;
                        IN.minY = nPX{a}.minY; IN.maxY = nPX{a}.maxY; IN.revertflag = true;
                        try
                            nmW(nmW==0)=NaN; nmW = nk_PerfScaleObj(nmW', IN)'; nmW(isnan(nmW))=0;
                        catch
                            fprintf('problem');
                        end
                    end

                case {'reducedim','remvarcomp',}
                    
                    if isfield(nPX{a}.mpp,'vec')
                       redvec = nPX{a}.mpp.vec;
                    elseif isfield(nPX{a}.mpp,'factors')
                       redvec = nPX{a}.mpp.factors{1};
                    end
                    
                    if isfield(nPX{a},'ind0')
                        ind0 = nPX{a}.ind0;
                        DR = nPX{a}.DR;
                    else
                        ind0 = 1:size(redvec,2);
                        DR = nPREPROC.ACTPARAM{a}.DR;
                    end
                    
                    mpp.vec = redvec(:,ind0);
                    
                    switch DR.RedMode                                        
                        case {'PCA','RobPCA','SparsePCA'}
                            if strcmp(DR.RedMode,'RobPCA'), 
                                DRsoft = 1; 
                            else
                                DRsoft = DR.DRsoft;
                            end
                            % Project back to input space
                            switch DRsoft
                                case 0
                                    nmW = reconstruct_data(nmW, mpp);
                                    nmP = logical(reconstruct_data(nmP, mpp));
                                case 1
                                    nmW = mpp.vec * nmW;
                                    nmP = logical(mpp.vec * nmP);
                            end
                        case 'NNMF'
                             nmW = mpp.vec * nmW;
                             nmP = logical(mpp.vec * nmP);
                        otherwise
                             nmW = reconstruct_data(nmW, nPX{a}.mpp);
                             nmP = logical(reconstruct_data(nmP, nPX{a}.mpp));
                    end

                    % If features had been removed prior to dimensionality
                    % reduction take care that you recover the original space
                    if isfield(nPX{a},'indNonRem') && ~isempty(nPX{a}.indNonRem) && sum(~nPX{a}.indNonRem) > 0
                        tmW = zeros(size(nPX{a}.indNonRem')); tmP = zeros(size(nPX{a}.indNonRem'));
                        tmW(nPX{a}.indNonRem) = nmW; nmW = tmW; tmP(nPX{a}.indNonRem) = nmP; nmP = tmP;  
                        clear tmW tmP;
                    end
                    reducedimfl = true;

                case {'elimzero','extfeat','extdim'}
                    %flg = true;
                    if isfield(nPX{a},'NonPruneVec')
                        IND = 'NonPruneVec';
                    elseif isfield(nPX{a},'indNonRem')
                        IND = 'indNonRem';
                    else
                        IND = 'ind';
                    end
                    try
                        tmW = zeros(numel(nPX{a}.(IND)),1); tmW(nPX{a}.(IND)) = nmW; nmW = tmW;   
                        tmP = false(numel(nPX{a}.(IND)),1); tmP(nPX{a}.(IND)) = nmP; nmP = tmP;
                        if ~isempty(PA)
                            tmPA = zeros(numel(nPX{a}.(IND)),1);
                            tmPA(nPX{a}.(IND)) = nmPA; nmPA = tmPA;
                        end
                    catch
                        fprintf('problem')
                    end
            end

        end
    
        if size(nmW,2) > 1, nmW = nmW'; nmP = nmP'; 
            if ~isempty(PA), nmPA = nmPA'; end
        end
        
        % Run univariate correlation analysis if no factorization method
        % had been applied to the data
        if ~decompfl(n)
            tY = zeros(size(Y,1),numel(nmP)); 
            tY(:,nmP) = Y(:,VI(Fu) == n); 
            % Rows with Nans will be removed by nk_CorrMat
            nmR  =   nk_CorrMat(L, tY,'pearson'); 
            nmSR =   nk_CorrMat(L, tY,'spearman');
        else
            nmR  =  nan(size(nmW)); nmSR = nmR;
        end
        
        if nM == 1 && mM>1
            for m=1:numel(varind)
                if procfl
                    indX = inp.X.dimvecx(m)+1:inp.X.dimvecx(m+1);
                    mW{m} = nmW(indX);
                    mP{m} = nmP(indX);
                    mR{m} = nmR(indX);
                    mSR{m} = nmSR(indX);
                    if ~isempty(PA), mPA{m} = nmPA(indX); end
                else
                    mW{m} = nmW;
                    mP{m} = nmP;
                    mR{m} = nmR;
                    mSR{m} = nmSR;
                    if ~isempty(PA), mPA{m} = nmPA; end
                end
            end
        else
            mW{n} = nmW;
            mP{n} = nmP;
            mR{n} = nmR;
            mSR{n} = nmSR; 
            if ~isempty(PA), mPA{n} = nmPA; end
        end
    else
        % Concatenate weight vectors
        nmW = zeros(numel(fVI),1); nmW(nmP) = W(lFuVI);
        mW = [mW; nmW]; 
        if ~isempty(PA), 
             nmPA = zeros(numel(fVI),1); nmPA(nmP) = PA(lFuVI);
             mPA = [mPA; nmPA];
        end
    end
end

for n=1:numel(mW)
    ind0 = mW{n}==0;
    mW{n}(ind0) = NaN;
end
