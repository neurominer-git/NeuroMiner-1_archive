function [rs, ds] = nk_GetTestPerf_matLRN(~, tXtest, ~, md, ~, ~)
global MODEFL

p = md.predict(md, tXtest);
switch MODEFL
    case 'regression'
        rs = p; ds = p;
    case 'classification'
        if isstruct(p)
            rs = p.yhat; ds = p.D;
        else
            rs = p ; ds = p;
        end
end

