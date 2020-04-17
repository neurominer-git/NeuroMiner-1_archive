function Mx = display_classcomparison_matrix(pvals, ticklabels, mw, sw, T)
% =========================================================================
% function Mx = display_classcomparison_matrix(pvals, ticklabels, mw, sw)
% =========================================================================
% Displays P value matrix of pair-wise classifier comparisons 
% and mean(sd) values of classifiers' performance
% 
% Inputs:
% -------
% pvals :           P value matrix
% ticklabels :      Classifier descriptions
% mw :              mean CV2 performance of classifiers (vector)
% sw :              sd CV2 performance of classifiers (vector)
%
% Outputs:
% --------
% Mx :              discretized classifier matrix
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 03/2020

sz = get(0,'ScreenSize');
win_wdth = sz(3)/2; win_hght = sz(4); win_x = sz(3)/2 - win_wdth/2; win_y = sz(4)/2 - win_hght/2;
if ~exist('T','var') || isempty(T)
    T = -log10([0.05 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 realmin]) ;
end
if istable(pvals)
    pval_mat = pvals.Variables;
else
    pval_mat = pvals;
end
% Prepare P value matrix
pval_mat(pval_mat==0)=realmin; M = -log10(pval_mat);
Mx = M;
Mx ( M < T(1) ) = NaN;
for i = 1:numel(T)-1
    Mx ( M>= T(i) & M<T(i+1) ) = i;
end

figure('Name','Comparative Model Analysis', 'Position',[win_x win_y win_wdth win_hght]); ax = axes; imagesc(Mx, 'AlphaData', ~isnan(Mx));
ax.Position = [0.15 0.15 0.775 0.55];
ax.CLim=[ 1 numel(T) ];
CLabels = cellstr(num2str(T','%1.2f')); CLabels{end} = sprintf('>%1.2f',T(end-1));
colormap(jet(numel(T)-1));c = colorbar('Ticks', 1:numel(T), 'TickLabels', CLabels); c.Label.String = '-log10(P value)';

if exist('ticklabels','var') && ~isempty(ticklabels) && numel(ticklabels)==size(pvals,1)
    ax.XTick = 1: numel(ticklabels);
    ax.XTickLabel = ticklabels;
    ax.XTickLabelRotation = 45;
    ax.YTick = ax.XTick;
    ax.YTickLabel = ticklabels;
    ax.TickLabelInterpreter='none';
end

if exist('mw','var') && exist('sw','var')
    bx = axes('Position',[0.15 0.725 0.7 0.25]); 
    errorbar(bx, mw, sw, 'LineStyle', 'none', 'Color', 'k'); hold on
    plot(bx, mw, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'MarkerEdgeColor', 'w');
    xlim([0.5 numel(mw)+0.5]);
    bx.XTickLabel = [];
    ylabel('Performance');
end