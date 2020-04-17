function [h, n_mat] = nk_PlotCobWeb(mat, groupnames, ha)

cla(ha);
[M,N] = size(mat);
idx = eye(M,N);
n_mat = bsxfun(@rdivide, mat , sum(mat))*100;
mcl = n_mat(~idx);
Ax = cell(M,N); Ux = cell(M,N);
for i=1:M
    for j=1:N
        Ax{i,j} = sprintf('%s\n%s', groupnames{i},groupnames{j});
        Ux{i,j} = '';
    end
end
Ax = [Ax(~idx) Ux(~idx)];
% Plot cobweb graph
h = spider([ones(numel(mcl),1)*100/N mcl ],'Misclassification web',50*ones(numel(mcl),1),Ax,[], ha);