function [rs, ds] = nk_GetTestPerf_LIBSVM(~, tXtest, Ytest, md, ~, ~)
global SVM LIBSVMPREDICT

[err_test, ~, predict_test] = feval( LIBSVMPREDICT, Ytest, tXtest, md, sprintf(' -b %g',SVM.LIBSVM.Optimization.b)); 
rs = err_test; ds = predict_test(:,1);


