function Perf = nk_MultiPerfQuant(expected, predicted, modus)

nsubj = numel(expected);
ngroups = unique(expected);
nG = numel(ngroups);

switch modus
    case 0
        errs = predicted ~= expected;
        Perf= ( 1 - sum(errs)/ nsubj) * 100;
    case 1
        Perfi = zeros(nG,1); 
        for i = 1:nG
            expi=ones(size(expected)); predi=ones(size(expected));
            indp = predicted ~= i; inde = expected ~= i;
            expi(inde)=-1; predi(indp)=-1;
            Perfi(i) = BAC(expi,predi);
        end
        Perf = mean(Perfi);
end

end