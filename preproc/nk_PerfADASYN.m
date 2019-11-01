function [Y, L] = nk_PerfADASYN(Y, L, IN)

adasyn_beta                     = [];   %let ADASYN choose default
adasyn_kDensity                 = [];   %let ADASYN choose default
adasyn_kSMOTE                   = [];   %let ADASYN choose default
adasyn_normalized               = false;    %false lets ADASYN handle normalization

L(L==-1)=0;

if exist('IN','var') && ~isempty('IN')
    if isfield(IN,'beta') && ~isempty(IN.beta), adasyn_beta = IN.beta; end
    if isfield(IN,'kDensity') && ~isempty(IN.kDensity), adasyn_kDensity = IN.kDensity; end
    if isfield(IN,'kSMOTE') && ~isempty(IN.kSMOTE), adasyn_kSMOTE = IN.kSMOTE; end
    if isfield(IN,'normalized') && ~isempty(IN.normalized), adasyn_normalized = IN.normalized; end
end

[Ysyn, Lsyn] = ADASYN(Y, L, adasyn_beta, adasyn_kDensity, adasyn_kSMOTE, adasyn_normalized);

Y=[Y;Ysyn]; L=[L;Lsyn];
L(~L)=-1;