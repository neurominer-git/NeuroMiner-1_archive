% =========================================================================
% FORMAT param = nk_OOCV_config(param, res)
% =========================================================================
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07 / 2011

function [param, act] = nk_OOCV_config(param)
global NM

meanflag = 1;
preproc = 0;
savemodels = 2;
saveoocvdata = 1;
if isfield(param,'MULTI') && param.MULTI.flag, 
    multiflag = true; groupmode = 2; 
else
    multiflag = false;
    groupmode = 1; 
end
trainwithCV2Ts = 1;

if ~isfield(NM.TrainParam,'OOCV')
    NM.TrainParam.OOCV.meanflag = meanflag;
    NM.TrainParam.OOCV.preproc = preproc;
    %NM.TrainParam.OOCV.savemodels = savemodels;
    NM.TrainParam.OOCV.groupmode = groupmode;
    NM.TrainParam.OOCV.multiflag = multiflag;
    NM.TrainParam.OOCV.trainwithCV2Ts = trainwithCV2Ts;
    %NM.TrainParam.OOCV.saveoocvdata = saveoocvdata;
end
o = nk_GetParamDescription2(NM,NM.TrainParam,'oocv');

nk_PrintLogo

if multiflag && ~strcmp(NM,'regression')
    act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> OOCV settings',0,'mq', ...
                    ['Data calibration [ ' o.preprocstr ' ]|'  ...
                     'Ensemble mode [ ' o.meanflagstr ' ]|' ...
                     'Group processing mode [ ' o.groupmodestr ' ]|' ...
                     'Model retraining mode [ ' o.trainwithCV2Tsstr ' ]'],1:4);
                     %'Save retrained models [ ' o.savemodelstr ' ]' ...
                     %'Save preprocessed independent test data [ ' o.saveoocvdatastr ' ]'],1:6);
else
    act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> OOCV settings',0,'mq', ...
                    ['Preprocessing [ ' o.preprocstr ' ]|'  ...
                     'Ensemble mode [ ' o.meanflagstr ' ]|' ...
                     'Model retraining mode [ ' o.trainwithCV2Tsstr ' ]'],[1 2 4]);
                     %'Save retrained models [ ' o.savemodelstr ' ]|' ...
                     %'Save preprocessed independent test data [ ' o.saveoocvdatastr ' ]'],[1 2 4 5 6]);
end

switch act

    case 1
        preproc = nk_input('Adjust OOCV data for mean effects',0, 'm', ...
                                ['No adjustment|' ...
                                 'Mean adjustment using mean function|' ...
                                 'Mean adjustment using pinv'],0:2, preproc);
        if preproc == 2
            preprocind = nk_input('Define index vector(s)',0,'e',[], size(NM.OOCV{1}.Y{1},1));
        end
        v
    case 2
        meanflag = nk_input('Ensemble construction mode', 0, 'm', ...
                                ['Aggregate ALL base learners into ONE ensemble|' ...
                                 'Compute mean decisions of CV1 partion ensembles before aggregating'], ...
                                 [1,2], meanflag);
    case 3
        groupmode = nk_input('Define group processing mode',0, 'm', ...
                                ['OOCV prediction at binary predictors'' optimum parameters|' ...
                                 'OOCV prediction at multi-group predictors'' optimum parameters|' ...
                                 'OOCV prediction at binary & multi-group predictors'' optimum parameters'], ...
                                 1:3, groupmode);
            
    case 4
        trainwithCV2Ts = nk_input('Model retraining mode', 0, 'm', ...
                                 ['As defined for CV training phase|' ...
                                  'Use all available CV data (CV1+CV2)'], ...
                                  [1,2], trainwithCV2Ts);
    
%     case 5
%         savemodels = nk_input('Save retrained models', 0, 'yes|no',[1,2], savemodels);
%     case 6
%         saveoocvdata = nk_input('Save retrained models', 0, 'yes|no',[1,2], saveoocvdata);
  
end

param.OOCV.preproc = preproc;
if exist('preprocind','var')
    param.OOCV.preprocind = preprocind;
end
param.OOCV.meanflag         = meanflag;
param.OOCV.groupmode        = groupmode;
param.OOCV.trainwithCV2Ts   = trainwithCV2Ts;
param.OOCV.savemodels       = savemodels;
param.OOCV.saveoocvdata     = saveoocvdata;

end