function performance = nk_MultiClassAssessConfMatrix(confmatrix, label, pred, errs)

sumconf=sum(sum(confmatrix));
P=0; inan = isnan(pred); pred(inan) = []; label(inan)=[];
% Calculate sub confusion matrix for EACH GROUP vs the REST
uL = unique(label);
for i=1:size(confmatrix,1)
    labeli = label; labeli(label~=uL(i)) = 0;labeli(label==uL(i)) = 1;
    predi = pred; predi(pred~=uL(i)) = 0; predi(pred==uL(i)) = 1;
	h = confmatrix(i,:);
	h(i) = []; FN = sum(h);
	TP = confmatrix(i,i);
	FP = sum(confmatrix(:,i)) - TP;
	TN = sumconf-TP-FP-FN;
	performance.class{i}.TP     = TP;
	performance.class{i}.TN     = TN;
	performance.class{i}.FP     = FP;
	performance.class{i}.FN     = FN;
    performance.class{i}.acc    = (TP+TN)/(TP+TN+FP+FN)*100;
	performance.class{i}.sens   = (TP / (TP + FN)) * 100; sens = TP / (TP + FN);
	performance.class{i}.spec   = (TN / (TN + FP)) * 100; spec = TN / (TN + FP);
    performance.class{i}.FPR    = (FP / (TN + FP)) * 100;
	performance.class{i}.PPV    = (TP / (TP + FP)) * 100;
	performance.class{i}.NPV    = (TN / (TN + FN)) * 100;
    performance.class{i}.AUC    = fastAUC(labeli,predi);
    TPrate = TP / ( TP + FN); TNrate = TN / ( TN + FP); 
    performance.class{i}.GMean  = sqrt(TPrate + TNrate);
	performance.class{i}.MCC    = (TP*TN-FP*FN)/sqrt((TP+FN)*(FP+TN)*(TP+FP)*(FN+TN));	
    performance.class{i}.Fscore = (2 * TP) / (2 * TP + FP + FN) ;
    performance.class{i}.BAC    = (performance.class{i}.sens + performance.class{i}.spec) / 2;
    performance.class{i}.pLR    = performance.class{i}.sens / (100 - performance.class{i}.spec);
    performance.class{i}.nLR    = (100- performance.class{i}.sens) / performance.class{i}.spec;
    performance.class{i}.PSI    = performance.class{i}.PPV + performance.class{i}.NPV - 100;
    performance.class{i}.NNP    = 1 / (performance.class{i}.PSI / 100);
    performance.class{i}.NND    = 1 / (sens - (1 - spec ));
    performance.class{i}.Youden = ( sens + spec ) - 1;
    performance.class{i}.DOR    = sens/(1-spec)/((1-spec)/sens);
	P = P + TP;
end
performance.confusion_matrix    = confmatrix;
performance.errors              = errs;
performance.accuracy            = P*100/sumconf;
performance.BAC                 = 0;
for i = 1:numel(performance.class)
    performance.BAC = performance.BAC + performance.class{1}.BAC;
end
performance.BAC = performance.BAC / numel(performance.class);

return