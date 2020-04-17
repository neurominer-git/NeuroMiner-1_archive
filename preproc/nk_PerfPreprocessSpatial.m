function [sY, sYocv, sCocv, inp, optfl, ukbin, uBINMOD, BINMOD] = ...
    nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, paramfl, BINMOD, kbin, ukbin)

global PREPROC
sYocv = []; sCocv = [];
%Eventually, apply spatial operations on image data
if isfield(paramfl,'PREPROC') && isfield(paramfl.PREPROC,'SPATIAL') && paramfl.PREPROC.SPATIAL.cubetype>1
    % Filter the data (imaging only)
    sY = cell(kbin,1); 
    if ~isempty(Yocv), sYocv = cell(kbin,1); end
    if ~isempty(Cocv), sCocv = cell(kbin,1); end
    optfl = true; 
    ukbin = kbin; 
    uBINMOD = 0; BINMOD = 1;
    for u=1:kbin
        uPREPROC = nk_SetParamChain(paramfl, u, PREPROC);
        paramfl.P{u} = nk_ReturnParamChain(uPREPROC, 1); 
        fprintf('\nPredictor #%g: Smoothing Training & CV data',u)
        sY{u} = nk_PerfSpatFilt2( Y, uPREPROC, paramfl.PV ); 
        if isfield(inp.X,'Yw'), 
            fprintf('\nSmoothing weighting map')
            inp.Yw = nk_PerfSpatFilt2( inp.X.Yw, uPREPROC, paramfl.PV ); 
        end
        if ~isempty(Yocv), 
            fprintf('\nPredictor #%g: Smoothing independent test data (%s)',u, inp.desc_oocv);
            sYocv{u} = nk_PerfSpatFilt2( Yocv, uPREPROC, paramfl.PV ); 
        end
        if ~isempty(Cocv), 
            fprintf('\nPredictor #%g: Smoothing calibration data',u)
            sCocv{u} = nk_PerfSpatFilt2( Cocv, uPREPROC, paramfl.PV ); 
        end
    end   
else
    uBINMOD = BINMOD;
    optfl = false;
    sY = Y;
    if ~isempty(Yocv), sYocv = Yocv; end
    if ~isempty(Cocv), sCocv = Cocv; end
    if isfield(inp,'X') && isfield(inp.X,'Yw') && ~isfield(inp,'Yw'), inp.Yw = inp.X.Yw; end
end