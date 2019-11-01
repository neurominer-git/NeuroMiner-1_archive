function handles = display_comparator(handles)
 
if ~isfield(handles,'PerfCompWin')
    WinTag = 'PerfCompWin';
    h = findobj('Tag',WinTag);
    if isempty(h)
       handles = init_fig(handles, WinTag);
    else
       handles.PerfTab.Win = h; 
    end
end
h = handles.PerfTab.Win ;
set(0,'CurrentFigure',h)

function handles = init_fig(handles, WinTag)

handles.PerfTab.Win = figure('NumberTitle','off','Name','NM Performance Comparator', 'Tag' ,WinTag,'MenuBar','none');

cnt=0; analyses=[];

for i=1:numel(handles.NM.analysis)
    if handles.NM.analysis{i}.status
        cnt = cnt+1;
        %analyses{cnt} = handles.NM.analysis{i}.id;
        analysescnt{cnt} = num2str(cnt);
        data{cnt, 1} = handles.NM.analysis{i}.id;
        data{cnt, 2} = false;
    end
end

handles.PerfTab.analysisselect = uitable('units', 'normalized', 'position', [0.05, 0.125, 0.9, 0.7], ...
                                        'ColumnName', {'Analyses','Select'}, ... 
                                        'ColumnFormat',{'char','logical'},...
                                        'ColumnEditable', [false, true],...
                                        'ColumnWidth',{300, 'auto'}, ...
                                        'RowName', analyses,...
                                        'data', data);
handles.PerfTab.addedit     = uicontrol( 'Style','edit', ...
                                        'Units', 'normalized', ...
                                        'Position', [0.05, 0.85, 0.9, 0.06]);
handles.PerfTab.addedit_label = uicontrol('Style','text', ...
                                        'Units', 'normalized', ...
                                        'Position', [0.05, 0.92, 0.9, 0.04], ...
                                        'String','Enter external predictor outputs here:');
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
fcnt = 0; cnt=1; a = zeros(1,sum(AnalysisSelection)); nA=numel(a);
for i=1:numel(handles.NM.analysis)
    if handles.NM.analysis{i}.status, 
        fcnt=fcnt+1; 
    else
        continue; 
    end
    if AnalysisSelection(fcnt)
        a(cnt) = i;
        CV = handles.NM.analysis{i}.params.cv;
        if i==1
            [ix, jx] = size(CV.TrainInd);
            [iy, jy] = size(CV.cvin{1,1}.TrainInd);
            PARAM = handles.NM.analysis{i}.params.TrainParam.SVM.GridParam;
            [~,~,PARAMFUN] = nk_GetScaleYAxisLabel(handles.NM.analysis{i}.params.TrainParam.SVM);
        else
            if size(CV.TrainInd,1) ~= ix 
                errordlg('The cross-validation structures of the selected analyses have an unequal number of CV2 repetitions');
                return
            elseif size(CV.TrainInd,2) ~= jx
                errordlg('The cross-validation structures of the selected analyses have an unequal number of CV2 folds');
                return
            elseif PARAM ~= handles.NM.analysis{i}.params.TrainParam.SVM.GridParam
                errordlg('The optimization criteria must be equal across the selected analyses');
            end
        end
        cnt = cnt+1;
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

if any(nanO), recomp = true; else, recomp = false; end   

nA = sum(AnalysisSelection);
if mPadd>0,
    PNames = cellstr([repmat('ExtPred_',mPadd,1) num2str((1:mPadd)')])';
else
    PNames = [];
end

% Map external predictor to cross-validation structure 
% (implement mapping to binary dichotomizers in multi-class case)
for curclass=1:handles.nclass
    if mPadd>0
        G = zeros(ix*jx,mPadd+nA);
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
    else
        G = zeros(ix*jx,nA);
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
                AnalNames = [AnalNames {sprintf('%s-M%g', handles.NM.analysis{a(i)}.id, M(g))}];
            else
                AnalNames = [AnalNames {handles.NM.analysis{a(i)}.id}];
            end
            
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
                        ll=ll+1;
                    end
                    NodesCnt(:,1) = NodesCnt(:,2)+1;
                end
            else
                G(:,ig) = AnalG.bestTS{curclass}(:);
            end
            ig=ig+1;
        end
    end
    AnalNames = [PNames AnalNames];
    handles.comparator_stats{curclass}.PredictorNames = AnalNames;
    handles.comparator_stats{curclass}.PredictorPerformances = G;
    if numel(AnalNames)>2
        handles.comparator_stats{curclass} = quadetest(G, AnalNames, handles.PerfTab.fileseltext.String);
    else
        handles.comparator_stats{curclass} = wilcoxon(G(:,1), G(:,2), 0.05);
    end
    
end

%tbl.colnames = data_table(1,1:end);
%tbl.rownames = data_table(2:end,1);
%tbl.array    = data_table(2:end,2:end);
%[ERR, STATUS, fil, typ] = tbl2file(tbl, handles.PerfTab.fileseltext.String, 'Performance Table');

% if STATUS
%     msgbox(['Data successfully exported to file ' fil]);
% else
%     warndlg(sprintf('Data NOT successfully exported to file.\nError: %s', ERR.id));
% end
