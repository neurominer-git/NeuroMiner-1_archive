function paramfl = nk_PrepPreprocParams(PREPROC, paramfl, analysis, modalind, ll, curlabel)
global VERBOSE MULTI

paramfl.PREPROC = PREPROC;        
if analysis.Model.NumPreMLParams > 0
    
    % Retrieve optimised parameters from GDdims
    paramfl.PXfull = nk_ReturnParamChain(PREPROC, 1);
    % Search for parameter column of given modality
    indm = analysis.Model.ModalityVec{1} == modalind;
    % No parameter column found
    NumMLParam = analysis.Model.NumParamDims - analysis.Model.NumPreMLParams ;
    
    for u=1:numel(analysis.bestP)
             
        if ~any(indm)
            pnt_vec = NumMLParam(u) + 1;
        else
            pnt_vec = NumMLParam(u) + find(indm);
        end
        
        if MULTI.flag && paramfl.multiflag
            if iscell(analysis.multi_bestPpos)
                Ppos = analysis.multi_bestPpos{ll, curlabel};
            else
                Ppos = analysis.multi_bestPpos(ll, curlabel);
            end
            paramfl.PXopt{u} = analysis.Model.ParamCombs{u}(Ppos, pnt_vec);
        else
            if ~iscell(analysis.bestP{1})
                paramfl.PXopt{u} = analysis.bestP{u}(ll, pnt_vec, curlabel);
            else
                paramfl.PXopt{u} = analysis.bestP{u}{ll}(:, pnt_vec, curlabel);
            end
        end
        % Retrieve unique combinations
        paramfl.PXunique{u} = unique(paramfl.PXopt{u},'rows','stable');
    end
else
    if VERBOSE, 
        fprintf('\n'); cprintf('blue','No preprocessing parameter optimization included in trained model');
    end
end