% =========================================================================
% =                          ROC ANALYSIS                                 =
% =========================================================================
function [hroc, hroc_random] = display_roc(handles, targets, predictions, axeshdl, clafl, linewidth)

if ~exist('axeshdl','var') || isempty(axeshdl), axeshdl = 'axes2'; end 
if ~exist('clafl','var') || isempty(clafl), clafl = true; end 
if ~exist('linewidth','var') || isempty(linewidth), linewidth = 2; end 

%GraphType = get(handles.selYaxis,'Value');
h_class         = get(handles.popupmenu1,'Value');
h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
h_classlist     = get(handles.popupmenu1,'String');

axes(handles.(axeshdl)); 
if clafl, cla; end
hold on

cl = 'k';
for i = 1:size(targets,2)
    if strcmp(handles.modeflag,'classification') && strcmpi(h_classlist{h_class},'Multi-group classifier')
        if h_onevsall_val == 1
            cl = handles.colptin(i,:);
        else
            cl = handles.colptin(h_onevsall_val-1,:);
        end
    end   
    %indnan = ~isnan(predictions(:,i));
    [Xsvm, Ysvm] = perfcurve2(targets(:,i), predictions(:,i), 1);
   
    %rocout = roc2( [targets(:,i), predictions(:,i),1);
    hroc(i) = plot(handles.(axeshdl),Xsvm,Ysvm, 'Color', cl, 'LineWidth', linewidth); 
end
xlabel(handles.(axeshdl),'False positive rate'); ylabel(handles.(axeshdl),'True positive rate'); 
handles.axes2.XLabel.Color='k'; 
%handles.axes2.XLabel.FontSize=12;
%handles.axes2.YLabel.FontSize=12;
hroc_random = plot(handles.(axeshdl),[0 1],[0 1],'k-.');
