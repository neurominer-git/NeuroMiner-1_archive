% =========================================================================
% =                        CLASSIFICATION PLOTS                           =
% =========================================================================
function handles = display_classplot_oocv(h, handles)

% Prepare axes
axes(handles.axes1); cla(handles.axes1); hold on
handles.axes1.Position = handles.axes1pos_orig;
handles.axes38.Visible = 'off';  cla(handles.axes38);
MS = 6;
MSoocv = 10;

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

% Get index to OOCV data
oocvind = handles.selCVoocv.Value - 1;

% Check whether the labels are known
labels_known = handles.OOCVinfo.Analyses{handles.curranal}.labels_known(oocvind);

switch GraphType
    case {1,2,3}
        P_h         = handles.BinClass{h}.mean_predictions;
        P_oocv_h    = nm_nanmedian(handles.OOCV(oocvind).data.predictions{h},2);
    case 4
        % Majority voting probabilities
        P_h = handles.BinClass{h}.prob_predictions(:,1);
        P_oocv_h    = nm_nanmedian(handles.OOCV(oocvind).data.prob_predictions{h},2);
    case 5
        % Mean Manjority voting probabilities
        P_h = handles.BinClass{h}.CV2grid.mean_predictions;
        P_oocv_h    = nm_nanmedian(handles.OOCV(oocvind).data.prob_predictions{h},2);
    case 6
        P_h = handles.BinClass{h}.CV2grid.mean_predictions;
        P_oocv_h    = nm_nanmedian(handles.OOCV(oocvind).data.prob_predictions{h},2);
end

% Mark groups with color
if isfield(handles.BinClass{h},'ind2')
    id1 = handles.BinClass{h}.ind1; 
    id2 = handles.BinClass{h}.ind2; 
else
    id1 = handles.BinClass{h}.ind1; 
    id2 = ~handles.BinClass{h}.ind1; 
end
legvecn = false(1,5);
handlevecn = cell(1,5);
% Print CV data: Group 1
handlevecn{1} = dotdensity( 1 ,P_h(id1), ...
    'dotEdgeColor', handles.colptin{handles.BinClass{h}.groupind(1)}, ...
    'dotFaceColor',handles.colptin{handles.BinClass{h}.groupind(1)}, ...
    'dotSize',MS, ...
    'dotMarker','o', ...
    'medianLine', 'on'); 

legvecn(1)=1;

% Print CV data: Group 2
if numel(handles.BinClass{h}.groupind) == 2
    CLP = 'o'; CLR = handles.colptin{handles.BinClass{h}.groupind(2)};
else
    CLP = 'o'; CLR = 'k';
end

handlevecn{2} = dotdensity(3,P_h(id2), ...
    'dotEdgeColor', CLR, ...
    'dotFaceColor', CLR, ...
    'dotSize',MS, ...
    'dotMarker', CLP, ...
    'medianLine', 'on');

legvecn(2)=1;

if labels_known
    if isfield(handles.OOCV(oocvind).data,'MultiResults')
        label_oocv_h = handles.OOCV(oocvind).data.MultiResults{1}.BinLabels{h};
    else
        label_oocv_h = handles.OOCV(oocvind).data.BinResults{1}.BinLabels{h};
    end

    ind0 =  label_oocv_h ~=0 ; fid0_oocv = find(~label_oocv_h);
    tP_oocv_h = P_oocv_h(ind0);
    id1_oocv = label_oocv_h(ind0) == 1 & ~isnan(tP_oocv_h) ; fid1_oocv = find(label_oocv_h(ind0) == 1);
    id2_oocv = label_oocv_h(ind0) == -1 & ~isnan(tP_oocv_h) ; fid2_oocv = find(label_oocv_h(ind0) == -1);

     % Print independent sample prediction using dot density plots: Group 1
     if sum(id1_oocv)
        [handlevecn{3},~,N{1},X{1},Y{1}] = dotdensity( 2 ,tP_oocv_h(id1_oocv), ...
            'dotEdgeColor', handles.colptin{handles.BinClass{h}.groupind(1)}, ...
            'dotFaceColor',handles.colptin{handles.BinClass{h}.groupind(1)}, ...
            'dotSize',MSoocv, ...
            'dotMarker','o', ...
            'medianLine', 'on');
        legvecn(3) = 1;
        fN{1} = fid1_oocv(N{1});
     end
    
     % Print independent sample prediction using dot density plots: Group 2
    if sum(id2_oocv)
        if numel(handles.BinClass{h}.groupind) == 2
            CLP = 'o'; CLR = handles.colptin{handles.BinClass{h}.groupind(2)};
        else
            CLP = 'o'; CLR = 'k';
        end
        [handlevecn{4},~,N{2},X{2},Y{2}] = dotdensity( 4, tP_oocv_h(id2_oocv), ...
            'dotEdgeColor', CLR, ...
            'dotFaceColor', CLR, ...
            'dotSize',MSoocv, ...
            'dotMarker', CLP, ...
            'medianLine', 'on');
        legvecn(4)=1;
        fN{2} = fid2_oocv(N{2});
    end

    if sum(~ind0)
        [handlevecn{5},~,N{3},X{3},Y{3}] = dotdensity( 5 , P_oocv_h(~ind0), ...
            'dotEdgeColor', 'k', ...
            'dotFaceColor', 'k', ...
            'dotSize',MSoocv, ...
            'dotMarker', 'o', ...
            'medianLine', 'on');
        legvecn(5)=1;
        fN{3} = fid0_oocv(N{3});
    end
    N=cell2mat(fN');X=cell2mat(X');Y=cell2mat(Y');
else
    [handlevecn{5},~,N,X,Y] = dotdensity( 5 , P_oocv_h, ...
        'dotEdgeColor', 'k', ...
        'dotFaceColor', 'k', ...
        'dotSize',MSoocv, ...
        'dotMarker', 'o', ...
        'medianLine', 'on');
    legvecn(5)=1;
end
    
% Define textbox info data 
pss = cell(1,numel(N));psslen=0;
for i=1:numel(pss)
    expgroupi = 'not labeled';
    if labels_known
        if label_oocv_h(i) > 0, 
            expgroupi = handles.BinClass{h}.groupnames{1}; 
        elseif  label_oocv_h(i) < 0, 
            expgroupi = handles.BinClass{h}.groupnames{2}; 
        end
    end
    if sign(P_oocv_h(i)) > 0,
        predgroupi = handles.BinClass{h}.groupnames{1}; 
    else
        predgroupi = handles.BinClass{h}.groupnames{2}; 
    end
    pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g'], handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind}{i}, expgroupi, predgroupi, P_oocv_h(i));
    if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.cases       = handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind};
figdata.x           = X;
figdata.y           = Y;
figdata.patterntext = pss(N);
figdata.parentui    = handles.pnBinary;
figdata.pnpos       = handles.pnBinary.Position;
figdata.figpos      = handles.figure1.Position;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... %'VerticalAlign', 'Top', ...
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);
    
if GraphType >3 
    probfx = 0.5;
    P_hx = P_h-0.5; handles.BinClass{h}.P_h = P_hx;
else
    probfx = 0; handles.BinClass{h}.P_h = P_h;
end

LegVec{1} = sprintf('Discovery: %s', handles.BinClass{h}.groupnames{1});
LegVec{2} = sprintf('Discovery: %s', handles.BinClass{h}.groupnames{2});
LegVec{3} = sprintf('OOCV: %s', handles.BinClass{h}.groupnames{1});
LegVec{4} = sprintf('OOCV: %s', handles.BinClass{h}.groupnames{2});
LegVec{5} = sprintf('OOCV: unlabeled');
handlevec = handlevecn(legvecn==1);
Hvec =[]; for i = 1:numel(handlevec), Hvec = [Hvec handlevec{i}]; end

legendvec = LegVec(legvecn==1);
    
switch GraphType

    case {1,2,3}
        switch handles.params.TrainParam.SVM.prog
            case {'MikSVM','MKLRVM'}
                algostr = 'RVM probability';
            case 'LIBSVM'
                algostr = 'SVM score';
            case 'MVTRVR'
                algostr = 'RVR score';
            case 'MEXELM'
                algostr = 'ELM score';
            case 'LIBLIN'
                switch handles.params.TrainParam.SVM.LIBLIN.classifier
                    case {0,6}
                        algostr = 'LR';
                    otherwise
                        algostr = 'SVM';
                end
                switch handles.params.TrainParam.SVM.LIBLIN.b
                    case 0
                        algostr = [algostr ' Score'];
                    case 1
                        algostr = [algostr ' Probability'];
                end
            case 'matLRN'
                algostr = sprintf('matLearn [ %s ]',char(handles.params.TrainParam.SVM.matLRN.algo));
            otherwise
                algostr = [handles.params.TrainParam.SVM.prog 'score'];
        end
    case 4
        algostr = 'OOT-Probability';
    case 5
        algostr = 'Mean OOT-Probability (95%-CIs)';
    case 6
        algostr = 'Mean OOT-Probability (SD)';
end
hx(1) = xlabel('Samples'); 
set(hx(1), ...%'FontSize',handles.AxisLabelSize-2,
    'FontWeight',handles.AxisLabelWeight);
hx(2) = ylabel(['Binned ' algostr]); 
set(hx(2), ...%'FontSize',handles.AxisLabelSize-2, ...
    'FontWeight',handles.AxisLabelWeight);
handles.legend_classplot = legend(Hvec, legendvec, 'Location','Best', 'FontUnits','normalized', 'FontSize', 8,'LineWidth',1);%,'FontSize',handles.LegendFontSize); 
legend('boxon')
flg = 'off'; flg2='off';  xlims = 6;
switch labels_known
    case 1
        flg='on';
        %% Extract contigency structure
        if isfield(handles.OOCV(oocvind).data,'BinResults')
            contigmat = handles.OOCV(oocvind).data.BinResults{1}.contingency{h};
        else
            contigmat = handles.OOCV(oocvind).data.MultiResults{1}.contingency{h};
        end
        confmatrix = [[contigmat.TP contigmat.FN]; [contigmat.FP contigmat.TN]];

        %% Display Contigency data
        %if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
        handles = display_contigmat(handles, contigmat);

        %% Display contingency plot
        handles.h_contig = display_contigplot(handles, confmatrix);

        if sum(ind0)>2 %&& numel(unique(label_oocv_h(ind0)))>1 
            %% Display ROC
            if isfield(contigmat,'AUC')
                flg2 = 'on'; [handles.hroc, handles.hroc_random] = display_roc(handles, label_oocv_h , P_oocv_h);
            end
            %% Display pie charts
            [handles.h1pie, handles.h2pie] = display_piecharts(handles, [], contigmat, label_oocv_h );   
        else
            xlims = 5;
        end
end
xlim(handles.axes1, [0 xlims]);
ylim(handles.axes1, 'auto');
handles.axes1.XTickLabel = [];
handles.pnRocCmds.Visible = flg2;
handles.pnPieCmds.Visible = flg2;
handles.pnContigCmds.Visible = flg;
handles.cmdExportCobWeb.Visible = flg;
handles.cmdMetricExport.Visible = flg;
handles.cmdExportAxes20.Visible = flg;

% Plot binary class deviding line
xLimits = get(handles.axes1,'XLim'); xLimitsVec = xLimits(1):xLimits(2);
zeroline = ones(1,numel(xLimitsVec))*probfx;
plot(handles.axes1,xLimitsVec,zeroline,'k--','LineWidth',handles.ZeroLineWidth)
    
