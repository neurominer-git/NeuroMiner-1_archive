function NM = nk_DefineOOCVCovars_config(NM)
    
    if isempty(NM.oocvvec)
        NM.covars_oocv = nk_DefineCovars_config(NM.n_subjects_all_oocv, NM.covars);
    else
        NM.covars_oocv             = NM.covars(NM.oocvvec,:);
        NM.covars(NM.oocvvec,:)   = [];
    end
end