function [xV, AnalP] = nk_GetAlgoWeightVec(SVM, Y, L, MD, decompfl, errorflag)
global STACKING

AnalP = []; xV=[];
switch SVM.prog
    case {'LIBSVM','LIBLIN'}
        switch SVM.kernel.kernstr
            case{' -t 0',' -t 4',' -t 5','lin', 'linear', 'none'} 
                %%%%%%%%%%%%% USE WEIGHTS OF MODEL %%%%%%%%%%%%
                xV = nk_GetPrimalW(MD); % Get weight vector over feature space
                if ~decompfl(1)
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
    case 'SEQOPT'
        nF = numel(STACKING.sel_anal);
        xV = zeros(nF,1);
        xV(MD.AnalSeq) = (MD.examsfreq(MD.examsfreq>0)/100)';
        
    otherwise
        if errorflag
            error(['Vizualisation for ' SVM.prog ' not supported yet. Please consult the NM Manual for more information']);
        end
end