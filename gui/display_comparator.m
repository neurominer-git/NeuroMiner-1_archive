function handles = display_comparator(handles, act)

if ~exist('act','var') || isempty(act), act = 'stats'; end
if ~isfield(handles,'PerfCompWin')
    WinTag = 'PerfCompWin';
    h = findobj('Tag',WinTag);
    if isempty(h)
       handles = init_fig(handles, WinTag, act);
    else
       handles.PerfTab.Win = h; 
    end
end
h = handles.PerfTab.Win ;
set(0,'CurrentFigure',h)
% -------------------------------------------------------------------------
function handles = init_fig(handles, WinTag, ActionMode)

cnt=0; analyses=[];
                                    
switch ActionMode
    
    case 'stats'
        
        for i=1:numel(handles.NM.analysis)
            if handles.NM.analysis{i}.status
                cnt = cnt+1;
                %analyses{cnt} = handles.NM.analysis{i}.id;
                analysescnt{cnt} = num2str(cnt);
                data{cnt, 1} = handles.NM.analysis{i}.id;
                data{cnt, 2} = false;
            end
        end
        
        handles.PerfTab.Win = figure('NumberTitle','off','Name','NM Statistical Performance Comparator', 'Tag' ,WinTag,'MenuBar','none');
        
        multiflag = false; togmult='off'; if numel(unique(handles.NM.label(~isnan(handles.NM.label))))>2, multiflag = true; togmult='on';end
        handles.PerfTab.analysisselect = uitable('units', 'normalized', 'position', [0.05, 0.110, 0.9, 0.7], ...
                                        'ColumnName', {'Analyses','Select'}, ... 
                                        'ColumnFormat',{'char','logical'},...
                                        'ColumnEditable', [false, true],...
                                        'ColumnWidth',{400, 'auto'}, ...
                                        'RowName', analyses,...
                                        'data', data);

        handles.PerfTab.addedit     = uicontrol( 'Style','edit', ...
                                        'Units', 'normalized', ...
                                        'Position', [0.05, 0.85, 0.7, 0.06], ...
                                        'TooltipString', sprintf(['Enter additional predictors here that were obtained outside of your NM analysis.' ...
                                                                '\nMake sure that prediction vectors have the same number of cases\n' ...
                                                                'and predictions match NM outputs.']));
        handles.PerfTab.addedit_label = uicontrol('Style','text', ...
                                        'Units', 'normalized', ...
                                        'Position', [0.05, 0.92, 0.7, 0.04], ...
                                        'String',sprintf('Enter external predictor outputs here [ %g cases, any number of predictors ]', size(handles.NM.cases,1)));
        handles.PerfTab.multiflag = uicontrol('Style','togglebutton', ...
                                        'units','normalized',...
                                        'Position',[0.80 0.85 0.15 0.06], ...
                                        'String','Multi-group', ...
                                        'Enable', togmult, ...
                                        'TooltipString', 'toggle button to compare multi-group or binary classifiers');   
        handles.PerfTab.fileseltext = uicontrol('Style','edit', ...
                                        'units','normalized', ...
                                        'Position',[0.05 0.04 0.58 0.06]);
        handles.PerfTab.fileseldlg  = uicontrol('Style','pushbutton', ...
                                        'units','normalized',...
                                        'Position',[0.64 0.04 0.15 0.06], ...
                                        'String','Save as', ...
                                        'Callback', {@saveas,handles});
        handles.PerfTab.tabulate    = uicontrol('Style','pushbutton', ...
                                        'units','normalized', ...
                                        'Position',[0.80 0.04 0.15 0.06], ...
                                        'String','Compare', ...
                                        'FontWeight', 'bold',...
                                        'BackgroundColor', rgb('lightblue'), ...
                                        'Callback', {@compare_predictors, handles});
    case 'visual'
        
        for i=1:numel(handles.NM.analysis)
            if handles.NM.analysis{i}.status
                cnt = cnt+1;
                %analyses{cnt} = handles.NM.analysis{i}.id;
                analysescnt{cnt} = num2str(cnt);
                data{cnt, 1} = handles.NM.analysis{i}.id;
                data{cnt, 2} = '';
                data{cnt, 3} = false;
            end
        end
        
        handles.PerfTab.Win = figure('NumberTitle','off','Name','NM Visual Performance Comparator', 'Tag' ,WinTag,'MenuBar','none');
        rowperfs=[];
        handles.PerfTab.analysisselect = uitable('units', 'normalized', 'position', [0.05, 0.110, 0.9, 0.45], ...
                                        'ColumnName', {'Analyses','Alias','Select'}, ... 
                                        'ColumnFormat',{'char','char','logical'},...
                                        'ColumnEditable', [false, true, true],...
                                        'ColumnWidth',{300, 'auto', 'auto'}, ...
                                        'RowName', analyses,...
                                        'data', data);
        switch handles.modeflag
            case 'classification'
                perfs1 = {'Balanced Accuracy', ...
                            'Accuracy', ...
                            'Sensitivity', ...
                            'Specificity', ...
                            'Positive Predictive Value', ...
                            'Negative Predictive Value', ...
                            'False Positive Rate', ...
                            'Prognostic Summary Index', ...
                            'Positive Likelihood Ratio', ...
                            'Negative Likelihood Ratio', ...
                            'Diagnostic Odds Ratio', ...
                            'Number Needed to Predict', ...
                            'Number Needed to Diagnose', ...
                            'Youden-Index', ...
                            'Matthews Correlation Coefficient', ...
                            'F-Score', ...
                            'G-Mean'}';
                        
                perfs2 = {'BAC', 'acc', 'sens', 'spec', 'PPV', 'NPV', 'FPR', 'PSI', 'pLR', 'nLR', 'DOR', 'NNP', 'NND', 'Youden', 'MCC', 'Fscore', 'Gmean'}';
                perfs3 = repmat({false},numel(perfs2),1);
                perfs4 = repmat({false},numel(perfs2),1);
                perfs5 = {'Red', 'LawnGreen', 'Blue', 'Gold', 'Orange', 'DarkSalmon', 'Gray', 'BurlyWood', 'DarkGreen', 'Turquoise', 'Fuchsia', 'MediumBlue', 'DarkGray', 'DarkRed', 'LightSkyBlue', 'DarkBlue', 'MediumPurple'}';
                perfs = [perfs1 perfs2 perfs3 perfs4 perfs5];
            case 'regression'
        end
        colors = rgb('getcolors');
        handles.PerfTab.perfselect     = uitable('units', 'normalized', 'position', [0.05, 0.58, 0.9, 0.40], ...
                                        'ColumnName', {'Performance Measures','','Select','Separate Axis', 'Color'}, ... 
                                        'ColumnFormat',{'char','char','logical','logical', colors.name'},...
                                        'ColumnEditable', [true, false, true, true, true],...
                                        'ColumnWidth',{300, 0, 'auto', 'auto', 50}, ...
                                        'RowName', rowperfs,...
                                        'data', perfs);
        handles.PerfTab.visualize    = uicontrol('Style','pushbutton', ...
                                        'units','normalized', ...
                                        'Position',[0.80 0.04 0.15 0.06], ...
                                        'String','Visualize', ...
                                        'FontWeight', 'bold',...
                                        'BackgroundColor', rgb('lightgreen'), ...
                                        'Callback', {@visualize_performances, handles});
        handles.PerfTab.sort         = uicontrol('Style','togglebutton', ...
                                        'units','normalized', ...
                                        'Position',[0.65 0.04 0.15 0.06], ...
                                        'String','Sort', ...
                                        'FontWeight', 'normal',...
                                        'BackgroundColor', rgb('lightgrey'));
        handles.PerfTab.line         = uicontrol('Style','togglebutton', ...
                                        'units','normalized', ...
                                        'Position',[0.50 0.04 0.15 0.06], ...
                                        'String','Line plots', ...
                                        'FontWeight', 'normal',...
                                        'BackgroundColor', rgb('lightgrey'));
end
 
function handles = saveas(src, evt, handles)

if ispc 
    ext = '*.xlsx';
else
    ext = '*.csv';
end

[FileName,PathName] = uiputfile(ext,'Save performance table','PerfTable');
handles.PerfTab.fileseltext.String = fullfile(PathName, FileName);


function handles = compare_predictors(src, evt, handles)

% Check whether path can be created
pth = fileparts(handles.PerfTab.fileseltext.String);
if isempty(handles.PerfTab.fileseltext.String) || isempty(pth), errordlg('Provide a valid output path before tabulating the data.'); return; end

curlabel = handles.curlabel;
curclass = 1;

if ~isfolder(pth)
    [status, msg] = mkdir(pth);
    if status
        rmdir(pth)
    else
        errordlg(msg)
        return
    end
end

AnalysisSelection = cell2mat(handles.PerfTab.analysisselect.Data(:,2));

if ~any(AnalysisSelection)
    errordlg('You have to select at least one analysis from the list')
    return
end

% Analyse cross-validation structures and optimization criteria
a = zeros(1,sum(AnalysisSelection)); nA=numel(a); fd = find(AnalysisSelection);
CV = handles.NM.analysis{fd(1)}.params.cv;
for i=1:nA
    a(i) = fd(i);
    if i==1
        [ix, jx] = size(handles.NM.analysis{fd(i)}.params.cv.TrainInd);
        [iy, jy] = size(handles.NM.analysis{fd(i)}.params.cv.cvin{1,1}.TrainInd);
        PARAM = handles.NM.analysis{fd(i)}.params.TrainParam.SVM.GridParam;
        [~,~,PARAMFUN] = nk_GetScaleYAxisLabel(handles.NM.analysis{fd(i)}.params.TrainParam.SVM);
    else
        if size(handles.NM.analysis{fd(i)}.params.cv.TrainInd,1) ~= ix 
            errordlg('The cross-validation structures of the selected analyses have an unequal number of CV2 repetitions');
            return
        elseif size(handles.NM.analysis{fd(i)}.params.cv.TrainInd,2) ~= jx
            errordlg('The cross-validation structures of the selected analyses have an unequal number of CV2 folds');
            return
        elseif PARAM ~= handles.NM.analysis{fd(i)}.params.TrainParam.SVM.GridParam
            errordlg('The optimization criteria must be equal across the selected analyses');
        end
    end
end

if ~isempty(handles.PerfTab.addedit.String)
    AdditionalPredictors = evalin('base', handles.PerfTab.addedit.String);
    [nPadd, mPadd] = size(AdditionalPredictors);
    if nPadd ~= size(handles.NM.label,1),
        errordlg('External prediction matrices should contain the same number of cases as the label data in your NM structure');
    end
    nanPadd = sum(isnan(AdditionalPredictors),2)>0;
else
    nPadd = size(handles.NM.label,1); mPadd = 0;
    nanPadd = false(nPadd,1);
end

% Determine cases with missing labels or data across selected analyses to
% find a population shared by all predictors

nanL    = sum(isnan(handles.NM.label),2)>0;
nanAnal = false(nPadd, nA);

for i=1:nA
    M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
    tNanAnal = false(nPadd, numel(M));
    for j=1:numel(M)
        tNanAnal(:,j) = sum(isnan(handles.NM.Y{M(j)}),2) == size(handles.NM.Y{M(j)},2);
    end
    nanAnal(:,i) = any(tNanAnal,2);
end

nanO = any([nanPadd nanL nanAnal],2);
Lg = handles.NM.label(:,curlabel);
Lg(nanO)=NaN;

% Do we have to recompute the prediction performanc because of cases with
% NaN labels or complete NaN data?
if any(nanO), recomp = true; else, recomp = false; end   

nA = sum(AnalysisSelection);
if mPadd>0,
    PNames = cellstr([repmat('ExtPred_',mPadd,1) num2str((1:mPadd)')])';
else
    PNames = [];
end

switch handles.PerfTab.multiflag.Value

    case 1
        G    = zeros(ix*jx, mPadd+nA, handles.ngroups);
        %Create one-vs-rest labels
        Lgfd = zeros(size(Lg,1),handles.ngroups);
        for curclass=1:handles.ngroups
             ind1 = Lg == curclass;
             ind2 = ~ind1;
             Lgfd(ind1,curclass) = 1; Lgfd(ind2,curclass)=-1;
        end
        Lgfd(nanO,:)=NaN;
        Gnames = cell(ix*jx,handles.ngroups);
        Gnames_Multi = cell(ix*jx,1);
        for curclass=1:handles.ngroups
            
             if mPadd>0
                %Compute performances for external predictors
                for g=1:mPadd
                    ll=1;
                    for f=1:ix
                        for d=1:jx
                             TsInd = CV.TestInd{f,d};
                             G(ll,g,curclass) = PARAMFUN(Lgfd(TsInd), AdditionalPredictors(TsInd));
                             ll=ll+1;
                        end
                    end
                end
             end
             
             lx = size(handles.NM.label,1); ig = mPadd+1;
             % now either get CV2 grids straightaway or recompute grid
             AnalNames = []; 
        
             for i=1:nA
                 
                AggrFlag = handles.NM.analysis{a(i)}.params.TrainParam.RFE.CV2Class.EnsembleStrategy.AggregationLevel;
                M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
                nGDdims = numel(handles.NM.analysis{a(i)}.GDdims);
                
                for g=1:nGDdims

                    AnalG = handles.NM.analysis{a(i)}.GDdims{g};
                    
                    if nGDdims > 1
                        AnalName = sprintf('%s_G%g-vs-REST_M%g', handles.NM.analysis{a(i)}.id, curclass, M(g));
                    elseif handles.nclass > 1
                        AnalName = sprintf('%s_G%g-vs-REST', handles.NM.analysis{a(i)}.id, curclass);
                    else
                        AnalName = handles.NM.analysis{a(i)}.id;
                    end
                    AnalName = regexprep(AnalName,'-','_');
                    if length(AnalName) > namelengthmax
                        warning('Variable name to long! Removing any underscores');
                        tAnalName = regexprep(AnalName,'_','');
                        if length(tAnalName) > namelengthmax
                            tAnalName = inputdlg(['The variable name is too long (max ' num2str(namelengthmax) ' characters. Please make manual adjustments'],'Error',[], AnalName);
                        end
                        AnalName = tAnalName;
                    end
                    
                    AnalNames = [AnalNames {AnalName}];

                    if recomp
                        
                        ll=1; 
                        NodesCnt = [ones(lx,1) zeros(lx,1)];

                        for f=1:ix

                            for d=1:jx
                                
                                TsInd = CV.TestInd{f,d};
                                if iscell(handles.NM.analysis{a(i)}.GDdims{g}.multi_bestPpos)
                                    nNodes = numel(handles.NM.analysis{a(i)}.GDdims{g}.multi_bestPpos{ll});
                                else
                                    nNodes = 1;
                                end
                                if ~AggrFlag
                                    nPred = iy*jy;
                                else 
                                    nPred = nNodes*iy*jy;
                                end
                                try
                                    NodesCnt(TsInd,2) = NodesCnt(TsInd,2) + nPred;
                                    llNodesCnt = NodesCnt(TsInd,:); N = numel(TsInd);
                                    Pred = AnalG.multi_predictions(TsInd, curlabel); 
                                    Pred = arrayfun( @(j) nm_nansum(Pred{j}(llNodesCnt(j,1):llNodesCnt(j,2))==curclass)*100/numel(llNodesCnt(j,1):llNodesCnt(j,2)), 1:N )';
                                    indP = Pred > 50;
                                    Pred(indP)=1; Pred(~indP)=-1;
                                catch 
                                    fprintf('problem')
                                end
                                G(ll,ig,curclass) = PARAMFUN(Lgfd(TsInd,curclass), Pred);
                                if i==1 && g==1, 
                                    Gnames{ll,curclass} = sprintf('CV2: R%g_F%g_C%g', f,d,curclass); 
                                    if curclass==1
                                        Gnames_Multi{ll} = sprintf('CV2: R%g_F%g', f,d); 
                                    end
                                end
                                ll=ll+1;
                            end
                            NodesCnt(:,1) = NodesCnt(:,2)+1;
                        end
                    else
                        G(:,ig,curclass) = AnalG.bestTS{curclass}(:);
                    end
                    ig=ig+1;
                end
             end
        end
        AnalNames = [PNames AnalNames];
        [ pth ,nam , ext ] = fileparts(handles.PerfTab.fileseltext.String);
        % Perform stats for each one-vs-rest classifier
        for curclass=1:handles.ngroups
            Filename = fullfile(pth, [nam sprintf('_G%g-vs-REST', curclass) ext]);
            handles.comparator_stats{curclass}.PredictorNames = AnalNames;
            handles.comparator_stats{curclass}.PredictorPerformances = G(:,:,curclass);
            if numel(AnalNames)>2
                handles.comparator_stats{curclass} = quadetest(G(:,:,curclass), Gnames(:,curclass), AnalNames, Filename);
            else
                handles.comparator_stats{curclass} = wilcoxon(G(:,1,curclass), G(:,2,curclass), 0.05);
            end
        end
        % Now perform multigroup test
        Gm = nm_nanmean(G,3);
        Filename = fullfile(pth, [nam '_Multi_Group' ext]);
        handles.comparator_stats_multi.PredictorNames = AnalNames;
        handles.comparator_stats_multi.PredictorPerformances = Gm;
        if numel(AnalNames)>2
            handles.comparator_stats{curclass} = quadetest(Gm, Gnames_Multi, AnalNames, Filename);
        else
            handles.comparator_stats_multi = wilcoxon(Gm(:,1), Gm(:,2), 0.05);
        end

    case 0
        % Map external predictor to cross-validation structure 
        % (implement mapping to binary dichotomizers in multi-class case)
        for curclass=1:handles.nclass
            G = zeros(ix*jx,mPadd+nA);
            if mPadd>0
                for g=1:mPadd
                    ll=1;
                    for f=1:ix
                        for d=1:jx
                            switch handles.modeflag
                                case 'classification'
                                    TsInd = CV.TestInd{f,d}(CV.classnew{f,d}{curclass}.ind);
                                    Lgfd = zeros(size(Lg,1),1);
                                    if numel(CV.classnew{f,d}{curclass}.groups)>1
                                        ind1 = Lg==CV.classnew{f,d}{curclass}.groups(1);
                                        ind2 = Lg==CV.classnew{f,d}{curclass}.groups(2);
                                    else
                                        ind1 = Lg==CV.classnew{f,d}{curclass}.groups;
                                        ind2 = ~ind1;
                                    end 
                                    Lgfd(ind1) = 1; Lgfd(ind2)=-1;
                                case 'regression'
                                    TsInd = CV.TestInd{f,d};
                                    Lgfd = Lg;
                            end
                            Lgfd = Lgfd(TsInd);
                            G(ll,g) = PARAMFUN(Lgfd, AdditionalPredictors(TsInd));
                            ll=ll+1;
                        end
                    end
                end
            end
            lx = size(handles.NM.label,1); ig = mPadd+1;

            % now either get CV2 grids straightaway or recompute grid
            AnalNames = []; Gnames = cell(ix*jx,1);

            for i=1:nA

                AggrFlag = handles.NM.analysis{a(i)}.params.TrainParam.RFE.CV2Class.EnsembleStrategy.AggregationLevel;
                M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
                nGDdims = numel(handles.NM.analysis{a(i)}.GDdims);

                for g=1:nGDdims

                    AnalG = handles.NM.analysis{a(i)}.GDdims{g};
                    if nGDdims > 1
                        AnalName = sprintf('%s_Cl%g_M%g', handles.NM.analysis{a(i)}.id, curclass, M(g));
                    elseif handles.nclass > 1
                        AnalName = sprintf('%s_Cl%g', handles.NM.analysis{a(i)}.id, curclass);
                    else
                        AnalName = handles.NM.analysis{a(i)}.id;
                    end
                    AnalName = regexprep(AnalName,'-','_');
                    if length(AnalName) > namelengthmax
                        warning('Variable name to long! Removing any underscores');
                        tAnalName = regexprep(AnalName,'_','');
                        if length(tAnalName) > namelengthmax
                            tAnalName = inputdlg(['The variable name is too long (max ' num2str(namelengthmax) ' characters. Please make manual adjustments'],'Error',[], AnalName);
                        end
                        AnalName = tAnalName;
                    end
                    AnalNames = [AnalNames {AnalName}];

                    if recomp

                        ll=1; 
                        NodesCnt = [ones(lx,1) zeros(lx,1)];

                        for f=1:ix

                            for d=1:jx

                                if iscell(handles.NM.analysis{a(i)}.GDdims{g}.bestPpos{curclass})
                                    nNodes = numel(handles.NM.analysis{a(i)}.GDdims{g}.bestPpos{curclass}{ll});
                                else
                                    nNodes = numel(handles.NM.analysis{a(i)}.GDdims{g}.bestPpos{curclass}(ll));
                                end
                                if ~AggrFlag
                                    nPred = iy*jy;
                                else 
                                    nPred = nNodes*iy*jy;
                                end
                                switch handles.modeflag
                                    case 'classification'
                                        TsInd = CV.TestInd{f,d}(CV.classnew{f,d}{curclass}.ind);
                                        Lgfd = zeros(size(Lg,1),1);
                                        if numel(CV.classnew{f,d}{curclass}.groups)>1
                                            ind1 = Lg==CV.classnew{f,d}{curclass}.groups(1);
                                            ind2 = Lg==CV.classnew{f,d}{curclass}.groups(2);
                                        else
                                            ind1 = Lg==CV.classnew{f,d}{curclass}.groups;
                                            ind2 = ~ind1;
                                        end 
                                        Lgfd(ind1) = 1; Lgfd(ind2)=-1;
                                    case 'regression'
                                        TsInd = CV.TestInd{f,d};
                                        Lgfd = Lg;
                                end
                                Lgfd = Lgfd(TsInd); 
                                NodesCnt(TsInd,2) = NodesCnt(TsInd,2) + nPred;
                                llNodesCnt = NodesCnt(TsInd,:); N = numel(TsInd);
                                Pred = AnalG.predictions(TsInd, curclass, curlabel); 
                                try
                                    Pred = arrayfun( @(j) nm_nanmedian(Pred{j}(llNodesCnt(j,1):llNodesCnt(j,2))), 1:N )';
                                catch 
                                    fprintf('problem')
                                end
                                G(ll,ig) = PARAMFUN(Lgfd, Pred);
                                if i==1 && g==1, Gnames{ll} = sprintf('CV2: R%g_F%g', f,d); end
                                ll=ll+1;
                            end
                            NodesCnt(:,1) = NodesCnt(:,2)+1;
                        end
                    else
                        ll=1;
                        for f=1:ix
                            for d=1:jx
                                G(ll,ig) = AnalG.bestTS{curclass}(f,d);
                                Gnames{ll} = sprintf('CV2: R%g_F%g', f,d);
                                ll=ll+1;
                            end
                        end
                    end
                    ig=ig+1;
                end
            end
            AnalNames = [PNames AnalNames];
            if handles.nclass > 1
                [ pth ,nam , ext ] = fileparts(handles.PerfTab.fileseltext.String);
                Filename = fullfile(pth, [nam sprintf('_Cl%g', curclass) ext]);
            else
                Filename = handles.PerfTab.fileseltext.String;
            end
            handles.comparator_stats{curclass}.PredictorNames = AnalNames;
            handles.comparator_stats{curclass}.PredictorPerformances = G;
            if numel(AnalNames)>2
                handles.comparator_stats{curclass} = quadetest(G, Gnames, AnalNames, Filename);
            else
                handles.comparator_stats{curclass} = wilcoxon(G(:,1), G(:,2), 0.05, Gnames, AnalNames, Filename);
            end

        end
end

function handles = visualize_performances(src, evt, handles)

Temp   = {'BAC', 'acc',       'sens', 'spec', 'PPV',    'NPV',        'FPR',  'PSI',       'pLR',        'nLR', 'DOR', 'NNP', 'NND', 'Youden', 'MCC', 'Fscore', 'Gmean'}'

AnalysisSelection = cell2mat(handles.PerfTab.analysisselect.Data(:,3));
AnalysisStrings = handles.PerfTab.analysisselect.Data(:,1);
AnalysisAliasStrings = handles.PerfTab.analysisselect.Data(:,2);
I = strcmp(AnalysisAliasStrings,'');
AnalysisAliasStrings(I) = AnalysisStrings(I);

PerfSelection = cell2mat(handles.PerfTab.perfselect.Data(:,3));
PerfFullStrings = handles.PerfTab.perfselect.Data(:,1);
PerfStrings = handles.PerfTab.perfselect.Data(:,2);
PerfSeparate = cell2mat(handles.PerfTab.perfselect.Data(:,4));
PerfColors = handles.PerfTab.perfselect.Data(:,5);

nA = numel(AnalysisSelection);
nAs = sum(AnalysisSelection);
nP = numel(PerfSelection);
nPs = sum(PerfSelection);
nS = sum(PerfSeparate);

Px = zeros(nAs, nPs);
cnt_i=0; 
curclass = 1;
gddims = 1;
for i=1:nA
    if AnalysisSelection(i)
        cnt_i=cnt_i+1;
        cnt_j=0;
        for j=1:nP
            if PerfSelection(j)
                cnt_j=cnt_j+1;
                Px(cnt_i,cnt_j) = handles.NM.analysis{i}.GDdims{gddims}.BinClass{curclass}.prob_contigency.(PerfStrings{j});
            end
        end
    end
end

AnalysisAliasStringsSel = AnalysisAliasStrings(AnalysisSelection);
% Eventually sort data according to means of rows
if handles.PerfTab.Win.Children(2).Value
    if nPs>1
        mPx = mean(Px,2);
    else
        mPx = Px;
    end
    [~,sI] = sort(mPx,'ascend');
    Px = Px(sI,:);
    AnalysisAliasStringsSel = AnalysisAliasStringsSel(sI);    
end

figure; hold on
if ~handles.PerfTab.Win.Children(1).Value
    if nPs>1
        if nS>0
            yyaxis left; 
            h{1} = bar(Px(:,~PerfSeparate(PerfSelection)),'grouped');
            yyaxis right;
            h{2} = bar(Px(:,PerfSeparate(PerfSelection)),'grouped');
        else
            h = bar(Px,'grouped');
        end
    else
        h = bar(Px);
    end
    ax = gca;
else
    if nS>0
        indLC = find(PerfSelection & ~PerfSeparate);
        indL = find(~PerfSeparate(PerfSelection)); 
        indRC = find(PerfSelection & PerfSeparate);
        indR = find(PerfSeparate(PerfSelection)); 
        yyaxis left
        for i=1:numel(indL)
            h{1,i} = plot(1:size(Px,1),Px(:,indL(i)),'-o','Color', rgb(PerfColors{indLC(i)}), 'MarkerFaceColor', rgb(PerfColors{indLC(i)}), 'LineWidth',2,'MarkerSize',7);
        end
        yyaxis right
        for i=1:numel(indR)
            h{2,i} = plot(1:size(Px,1),Px(:,indR(i)),'-X','Color', rgb(PerfColors{indRC(i)}),'MarkerFaceColor', rgb(PerfColors{indRC(i)}),'LineWidth',2,'MarkerSize',7);
        end
    else
        indSel = find(PerfSelection);
        for i=1:numel(indSel)
            h{i} = plot(1:size(Px,1),Px(:,i),'-o','Color',rgb(PerfColors{indSel(i)}), 'MarkerFaceColor', rgb(PerfColors{indSel(i)}), 'LineWidth',2,'MarkerSize',7);
        end
    end
    ax = gca;
    ax.XTick=1:size(Px,1);
end

legend(PerfFullStrings(PerfSelection),'Location','best'); 
ax.XTickLabel = AnalysisAliasStringsSel;
ax.XTickLabelRotation = 45;
ax.XLabel.String = 'Analyses';
ax.XLabel.FontWeight = 'bold';
ax.YLabel.String = 'Performance Measure(s)';
ax.YLabel.FontWeight = 'bold';