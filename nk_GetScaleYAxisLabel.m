function [yl, ylb] = nk_GetScaleYAxisLabel(SVM)
global PREPROC NM

switch SVM.GridParam
    case 1
        yl = [0 100];   ylb = 'Accuracy [%]';
    case 2
        yl = [0 100];   ylb = 'Sensitivity [%]';
    case 3
        yl = [0 100];   ylb = 'Specificity [%]';
    case 4
        yl = [0 100];   ylb = 'False Positive Rate [%]';
    case 5
        yl = [0 100];   ylb = 'Positive Predictive Value [%]';
    case 6
        yl = [-1 1];    ylb = 'Matthews Correlation Coefficient';       
    case 7
        yl = [0 1];     ylb = 'AUC';
    case 9
        yl = [0 100];   ylb = 'Mean squared error';
    case 10
        yl = [0 100];   ylb = 'Squared correlation coefficient [% explained variance]';
    case 11
        yl = [0 100];   ylb = 'Normalized root of mean squared deviation [%]';
    case 12
        yl = [0 1];     ylb = 'Root of mean squared deviation';
    case 13
        yl = [0 1];     ylb = 'Gmean';
    case 14
        yl = [0 100];   ylb = 'Balanced Accuracy [%]';
    case 15
        yl = [0 100];   ylb = 'F1-Score';
    case 16
        yl = [-1 1];    ylb = 'correlation coefficient';
    case 17
        yl = [0 100];   ylb = 'Enhanced Balanced Accuracy [%]';
    case 18
        if iscell(PREPROC),
            iPREPROC = PREPROC{1}; else, iPREPROC = PREPROC;
        end
        if isfield(iPREPROC,'LABELMOD') && isfield(iPREPROC.LABELMOD,'TARGETSCALE') && iPREPROC.LABELMOD.TARGETSCALE
            yl = [0 1]; 
        else
            if ~isempty(NM), 
                yl = [0 nk_Range(NM.label) ];
            else
                 yl = [0 nk_Range(evalin('base','NM.label'))];
            end
        end
        ylb = 'Mean Average Error [Label Range]';
    case 19
         yl = [-100 100];   ylb = 'Prognostic Summary Index [%]';
    case 20
         yl = [0 10]; ylb = 'Number Needed to Predict'; 
end

return