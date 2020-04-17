function [IN, optD, optCrit] = nk_OptPredSeq(D, L, IN, C, nCutOff, Lims, Crit, Ddesc)
% =========================================================================
% FORMAT function R = nk_OptPredSeq(D, L, nCutOff, Crit)
% =========================================================================
% This function optimizes the sequential combination of predictive models
% by ranking cases according to predictive scores and testing different 
% percentile thresholds for passing ranked cases on to the next predictive
% model
% Inputs:
% -------
% D :           Predictive score matrix, where each column contains the
%               predictive output of a model and each row the output of all 
%               models for a single case
% L :           The label vector
% IN :          A previously learned sequence predictor that should be
%               applied to new data
% C :           Combinations of predictive models to be tested. Each column
%               in C is a predictive model output stored in D 
% nCutOff [opt]:No. of thresholds to be tested x 2 (starting from 50%, 
% Lims [opt]:   lower and upper percentile limits
%               growing in both directions to Lims(2) / Lims(1)
% Crit [opt] :  an optimization criterion such as BAC
% Ddesc [opt] : a description of the nodes in the prediction sequence
%
% Output:
% -------
% IN :           Optimization results structure consisting of the following
%                fields:
%  IN.OPT :
%  IN.D :
%  IN.optD :
%  IN.Nremain :
%  IN.Crit :
%  IN.(Crit) :
%  IN.examsfreq :
%  IN.fcnt :
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 10/2018
global VERBOSE SVM

if isempty(IN)
    
    % Train sequence predictor
    if ~exist('nCutOff','var') || isempty(nCutOff), nCutOff = 5; end
    if ~exist('Crit','var') || isempty(Crit), Crit = 'BAC'; end
    if ~exist('Ddesc','var') || isempty(Ddesc)
        Ddesc = cellstr([repmat('Model ',size(D,2),1) num2str((1:size(D,2))')]); 
    elseif numel(Ddesc) ~= size(D,2)
        error('Model descriptor entries should correspond to the number of models in the model matrix');
    end

    nC = size(C,1);
    m = size(D,1);

    % Loop through a candidate prognostic workflows
    optCrit = []; cnt=1;
    for j=1:nC

        % Get model sequence 
        jC = C(j,~isnan(C(j,:)));

        if VERBOSE, fprintf('\n\n');cprintf('black*','Working on model sequence %g/%g:', j, nC); cprintf('blue','\t%s ', strjoin(Ddesc(jC),', ')); end

        % Get model outputs for given sequence in C
        jD = D(:,jC); nD = size(jD,2);
        
        % Build threshold vectors
        switch SVM.SEQOPT.AnchorType
            case 1
                Anchor = sum(jD < 0) * 100 / m; LL=zeros(numel(Anchor),1); UL=LL;
                for z=1:numel(Anchor)
                    Lthr = jD( jD(:,z) < 0, z );
                    pLthr = percentile(Lthr, 100 - Lims(1)); LL(z) = sum(jD(:,z) <= pLthr)*100/m;
                    Uthr = jD( jD(:,z) > 0, z ); 
                    pUthr = percentile(Uthr, Lims(2)); UL(z) = sum(jD(:,z) <= pUthr)*100/m;
                end
            case 2       
                Anchor = repmat(50,1,nD);
                LL = Anchor - Lims(1);
                UL = Anchor + Lims(2);
        end
        vecneg = cell(nD,1);
        vecpos = cell(nD,1);
        for z=1:nD
            vecneg{z} = Anchor(z): -((Anchor(z)-LL(z)) / nCutOff) : LL(z); vecneg{z}(1)=[];
            vecpos{z} = Anchor(z):  (UL(z)-Anchor(z)) / nCutOff : UL(z); vecpos{z}(1)=[];
        end
        
        % Model performance of first model in sequence
        OPT = ALLPARAM(L,jD(:,1));

        % Prepare results container for optimization
        tIN = struct('Sequence', jC, ...
                'OPT', OPT, ... 
                'fcnt', 1, ...
                'D', jD, ...
                'optD', jD(:,1), ...
                'Nremain', ones(m,1), ...
                'Crit', Crit);
        tIN.(Crit)(1) = OPT.(Crit);
        tIN.vecneg = vecneg;
        tIN.vecpos = vecpos;
        
        % Run sequential optimization
        for i=1:nD-1, 
            tIN = OptPredSeq(tIN, i, L, Crit); 
        end

        % Analyze frequencues of models selected in optimized workflow
        tIN.uN = unique(tIN.Nremain)';
        tIN.AnalSeq = jC(tIN.uN);
        tIN.AnalSeqDesc = Ddesc(jC(tIN.uN));
        tIN.examsfreq = zeros(1,numel(tIN.Sequence));
        for i=1:numel(tIN.Sequence)
            tIN.examsfreq(i) = sum(tIN.Nremain == i)*100/m;
        end
        if VERBOSE, fprintf('\nOptimized sequence performance: %1.2f', tIN.(Crit)(end)); end
        if j==1, optCrit = tIN.OPT.(Crit); IN = tIN; end

        if tIN.(Crit)(end) > optCrit(end)
            if VERBOSE, fprintf('\n');cprintf('blue*','+++ New optimal sequence: %s ', strjoin(Ddesc(jC),', ')); end
            IN = tIN; 
            cnt = cnt +1;
            optCrit(cnt) = IN.(Crit)(end);
        end
    end
    optD = IN.optD;
else
    % Apply trained sequence predictor to data
    D = D(:,IN.AnalSeq);
    Nremain = repmat(IN.uN(1),size(D,1),1);
    optD = D(:,1);
    for j=1:numel(IN.AnalSeq)-1
        % Apply absolute thresholds (questionable whether this is the best
        % strategy or alternatively apply learned percentiles)
        fI = find(Nremain==IN.uN(j));
        lthr = IN.optlthr(j); uthr = IN.optuthr(j);
        ind = D(fI,IN.uN(j)) >= lthr & D(fI,IN.uN(j)) <= uthr;
        Nremain(fI(ind)) = IN.uN(j+1);
        optD(fI(ind)) = D(fI(ind),j+1);
    end
    if exist('L','var') && ~isempty(L)
        OPT = ALLPARAM(L,optD);
        optCrit = OPT.(IN.Crit);
    end
end

function R = OptPredSeq(R, I, L, Crit)
global VERBOSE SVM

if isempty(VERBOSE), VERBOSE = true; end

optuvec = R.vecpos{I}(end);
optlvec = R.vecneg{I}(1);
optD = R.optD;
optRemain = R.Nremain;
fI = find(R.Nremain==I);
if isempty(fI) 
    fI = find(R.Nremain==I-1);
end
rOPT = R.OPT;
if VERBOSE, fprintf('\nProcessing: %g cases in predictive workflow node %g: ', numel(fI), I); end
DI = R.D(:,I);
warning off

for j=1:numel(R.vecneg{I})
    
    jD = R.optD;
    jDI = DI(fI); jDI(isnan(jDI))=[];
    lthr = percentile(jDI, R.vecneg{I}(j));
    uthr = percentile(jDI, R.vecpos{I}(j));
    if ~exist('optuthr','var')
         optlthr = lthr; optuthr = uthr;
    end
    ind = jDI >= lthr & jDI <= uthr;
    if ~any(ind), continue; end
    
    jNremain = R.Nremain;
    jNremain(fI(ind)) = I+1; 
   
    % Performance of new model in given sample defined by fI(ind)
    switch SVM.SEQOPT.Mode
        case 1
            prevOPT = ALLPARAM(L(fI(ind)),jD(fI(ind)));
        case 2
            prevOPT = rOPT;
    end
    jD(fI(ind)) = R.D(fI(ind),I+1);
    ijOPT = ALLPARAM(L, jD);
    try
        if ijOPT.(Crit) > prevOPT.(Crit)
            optlvec     = R.vecneg{I}(j);
            optuvec     = R.vecpos{I}(j);
            optlthr     = lthr; 
            optuthr     = uthr;
            optRemain   = jNremain;
            rOPT        = ijOPT;
            optD        = jD;
            if VERBOSE, fprintf('\nNext Node %g: %1.2f [%4g cases; Lower thresh: %1.2f, Upper thresh: %1.2f]', ...
                    I+1, rOPT.(Crit), numel(fI(ind)), lthr, uthr); end
        else
            if VERBOSE, fprintf('.'); end
        end
    catch
        fprintf('problem')
    end
end
R.optlvec(I) = optlvec;
R.optuvec(I) = optuvec;
R.optlthr(I) = optlthr;
R.optuthr(I) = optuthr;
R.optD = optD;
R.Nremain = optRemain;
R.jOPT(I) = rOPT;
R.(Crit)(I+1)= rOPT.(Crit);
R.OPT = rOPT;
