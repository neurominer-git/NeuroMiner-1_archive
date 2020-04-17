function p_map = compute_analytical_pvals(labels, Y, W)

%w = MODEL_SVM.model.SVs' * MODEL_SVM.model.sv_coef;
p       = sum(labels==1)/max(size(labels));
w       = W;
X       = Y;
r       = size(Y,1);
J       = ones(r,1);
K       = X*X';
Z       = inv(K)+(inv(K)*J*inv(-J'*inv(K)*J)*J'*inv(K));
C       = X'*Z;
%SD=sqrt(sum(C.^2,2));
SD      = sqrt(sum(C.^2,2)*(4*p-4*p^2));
mean_w  = sum(C,2)*(2*p-1);
p_map   = 2*normcdf(-abs(w-mean_w),zeros(size(w)),SD);

end