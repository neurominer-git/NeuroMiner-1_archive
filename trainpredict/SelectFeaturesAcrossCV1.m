function [ IN, OUT ] = SelectFeaturesAcrossCV1( IN, OUT, Param, OptMode )
global VERBOSE

nc = size(OUT.F,3);
ns = size(IN.F,3);

if check_fsizes(IN);
    cprintf('err','\nFeature masks have unequal dimensionalities due to feature preprocessing setup. Skipping Cross-CV1 feature operations.')
    return
end

for curclass=1:IN.nclass
    
    if VERBOSE, fprintf('\nAggregating feature masks for predictor #%g',curclass); end
    kIndMat = [];
    
    for i=1:IN.nperms % loop through CV1 permutations

        for j=1:IN.nfolds  % loop through CV1 folds
        
            %% Aggregate feature subspace mask :
            if nc > 1 
                % Every dichotomizer has its own feature subspace mask
                kIndMat = [kIndMat IN.F{i,j,curclass}(:,OUT.kxVec{i,j,curclass}(OUT.F{i,j,curclass}))];
            elseif nc == 1 && ns ==1
                kIndMat = [kIndMat IN.F{i,j}];
            else
                % All dichotomizers share one feature subspace mask
                % This mode is used for multi-group optimization
                kIndMat = [kIndMat IN.F{i,j}(:,OUT.kxVec{i,j,curclass}(OUT.F{i,j}))];  
            end
        end
    end
    
    % Probabilistic Feature Extraction across CV1 partitions
    tkIndMat = ProbabilisticFea(kIndMat,Param.PFE);
     
    % Overwrite existing feature masks with probabilistic mask
    ll=1;
    for i=1:IN.nperms % loop through CV1 permutations

        for j=1:IN.nfolds  % loop through CV1 folds
            
            if size(tkIndMat,2) > 1, kIndMat_ll = tkIndMat(:,ll); else, kIndMat_ll = tkIndMat; end
            if ~any(kIndMat_ll),
                error('Your feature selection procedure returned an empty feature space. Relax your selection parameters!'); 
            end
            if nc == 1 
                if ns > 1
                    IN.F{i,j,curclass} = kIndMat_ll;
                else
                    IN.F{i,j} = kIndMat_ll;
                end
            else
                OUT.F{i,j,curclass} = 1;
                OUT.Weights{i,j,curclass} = 1;
                IN.F{i,j,curclass} = kIndMat_ll;
            end    
            ll=ll+1;
        end
    end
end

% Retrain predictors using new feature masks
[IN, OUT] = FoldPerm(IN, OUT, 'Retrain for PFC across CV1 partitions', OptMode, 0, 0, Param.SubSpaceStepping); 
for curclass=1:IN.nclass

    % Overwrite existing feature masks with probabilistic mask
    for i=1:IN.nperms % loop through CV1 permutations

        for j=1:IN.nfolds  % loop through CV1 folds

            OUT.TrHDperf(i,j,curclass) = OUT.tr{i,j,curclass};
            OUT.TrHTperf(i,j,curclass) = OUT.tr{i,j,curclass};       
            OUT.CVHDperf(i,j,curclass) = OUT.ts{i,j,curclass};
            OUT.CVHTperf(i,j,curclass) = OUT.ts{i,j,curclass}; 

        end
    end
end

% _________________________________________________________________________
% Helper function that checks whether dimensionalities of training data 
% partitions are identical or not
function [uneqfl, u_dims, max_u_dims] = check_fsizes(IN)

n_dims = zeros(IN.nperms * IN.nfolds,1);
cnt = 1;
for i=1:IN.nperms % loop through CV1 permutations
   for j=1:IN.nfolds  % loop through CV1 folds
       n_dims(cnt) = numel(IN.F{i,j});
       cnt=cnt+1;
   end
end

u_dims = unique(n_dims);

if numel(u_dims)>1, 
    uneqfl = true; 
    max_u_dims = max(u_dims);
else
    uneqfl = false;
    max_u_dims = u_dims;
end

