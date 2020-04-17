function GDanalysis = nk_MultiPerfComp(GDanalysis, multi_pred, label, ngroups)

% Compute multi-class probabilties
multi_prob = nk_ConvProbabilities(multi_pred, ngroups);

lx = size(multi_prob,1);
pred = zeros(lx,1); 
stdpred = pred; ci1 = pred; ci2 = pred;
%numpred = zeros(lx,1);

% Loop through cases
if iscell(multi_pred)
    for i=1:lx
        if isempty(multi_pred{i})
            pred(i) = NaN; stdpred(i) = NaN; ci1(i) = NaN; ci2(i) = NaN; errs(i) = NaN;
        else
            %numpred(i) = length(multi_pred{i});
            % Maximum probability decides about multi-class membership
            [~,pred(i)] = max(multi_prob(i,:));
            % Is this useful: ?
            stdpred(i) = std(multi_pred{i});
            ci = percentile(multi_pred{i},[2.5 97.5]);
            ci1(i) = ci(1); ci2(i) = ci(2);
        end
    end
else
    %numpred = size(multi_pred,2);
    [~, pred] = max(multi_prob,[],2);
    stdpred = std(multi_pred,[],2);
    ci = cell2mat(arrayfun( @(i) percentile(multi_pred(i,:),[2.5 97.5]),1:lx,'UniformOutput',false )');
    ci1 = ci(:,1); ci2 = ci(:,2);
end

ind = ~isnan(pred);
if ~isempty(label)
    errs(ind) = label(ind)~= pred(ind);
    confmatrix = zeros(ngroups);
    % Compute confusion matrix
    for i=1:lx
        if ~ind(i), continue, end;
        confmatrix(label(i),pred(i)) = confmatrix(label(i),pred(i)) + 1;
    end
    % Compute performance measures and assign data to output
    GDanalysis.MultiClass = nk_MultiClassAssessConfMatrix(confmatrix, label, pred, errs);
end
GDanalysis.multi_probabilitiesCV2 = multi_prob;
GDanalysis.multi_predictionsCV2 = pred;
GDanalysis.multi_predictionsCV2_std = stdpred;
GDanalysis.multi_predictionsCV2_ci1= ci1;
GDanalysis.multi_predictionsCV2_ci2= ci2;

end