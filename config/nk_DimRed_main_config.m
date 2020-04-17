function [DR, PX, act] = nk_DimRed_main_config(DR, PX, parentstr, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
CALIBUSE = 2; 

if ~defaultsfl
    if isempty(DR), [DR, PX] = nk_DimRed_main_config(DR, PX, parentstr, true); end
    if isfield(DR,'RedMode'),   RedMode = DR.RedMode; else RedMode = []; end
    if isfield(DR,'dims'),      dims = DR.dims; else dims = []; end
    if isfield(DR,'PercMode'),  PercMode = DR.PercMode; else PercMode = []; end
    if isfield(DR,'CALIBUSE'),  CALIBUSE = DR.CALIBUSE; end
    if ~ischar(RedMode),
        DRSTR_REDMODE = 'undefined';
    else
        DRSTR_REDMODE = RedMode;
    end
    
    if strcmp(DRSTR_REDMODE,'PCA')
        if isempty(PercMode)
            DRSTR_PERCMODE = '(undefined)';
        else
            switch PercMode
                case 1
                    DRSTR_PERCMODE = ' (Absolute)';
                case 2
                    DRSTR_PERCMODE = ' (Percentage)';
                case 3
                    DRSTR_PERCMODE = ' (Energy)';
            end
        end
    else
        DRSTR_PERCMODE = '';
    end
    
    if isempty(dims)
        DRSTR_DIM = 'undefined';
    else
        DRSTR_DIM = [ nk_ConcatParamstr(dims) DRSTR_PERCMODE ];
    end
    if strcmp(DRSTR_REDMODE,'PLS')
        menustr = ['Select dimensionality reduction method [ ' DRSTR_REDMODE ' ]'];
        menuact = 1;
    else
        menustr = ['Select dimensionality reduction method [ ' DRSTR_REDMODE ' ]|', ...
                   'Define dimensionality of mapping results [ ' DRSTR_DIM ' ]'];
        menuact = 1:2;
    end
    
    [menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);
    
    nk_PrintLogo
    mestr = 'Dimensionality reduction'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    switch act
        case 1
            [DR, PX] = return_dimred(DR, PX, parentstr);
            [DR.dims, DR.PercMode] = nk_ExtDim_config(DR.RedMode, [], [], true);
        case 2 % for compatibility
            t_act = 1; while t_act > 0, [dims, PercMode, t_act] = nk_ExtDim_config(RedMode, PercMode, dims, 0, navistr);end
            DR.dims = dims; DR.PercMode = PercMode;
        case 1000
            CALIBUSE = nk_AskCalibUse_config(mestr, CALIBUSE); 
    end
else
    [DR, PX] = return_dimred(DR, PX, parentstr, true);
    [DR.dims, DR.PercMode] = nk_ExtDim_config(DR.RedMode, [], [], true);
    act = 0;
end

% Add parameters for dimensionality reduction optimization
PX = nk_AddParam(DR.dims, 'dimensions', 1, PX);
if numel(DR.dims)<2, DR.EXTOPT = true; else, DR.EXTOPT = false; end
DR.CALIBUSE = CALIBUSE;

end

function [DR, PX] = return_dimred(DR, PX, parentstr, defaultsfl)
global EXPERT NM
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl

   
    menustr = ['Principal Component Analysis                   (PCA)|' ...
               'Robust Principal Component Analysis            (significantly slower than PCA but more accurate in low N problems)|' ...
               'Non-negative Matrix Factorization (NNMF)       (significantly slower than PCA but more parsimonious)|' ...
               'Orthogonal Projective NNMF                     (significantly slower but more robust compared to NNMF)|' ... 
               'Fast-than-fast NNMF                            (fastest NNMF through random prejections and Nesterov iterations)|' ...
               'Partial Least Squares                          (PLS)|' ...
               'Sparse PCA                                     (... a ''statistical segmentation'' tool)|' ...
               'Probabilistic PCA                              (low-D only, >32 GB RAM for high-D data)|' ...
               'Deep Autoencoder                               (low-D only)|' ...
               'Factor analysis                                ()|' ...
               'Locality Preserving Projections                (low-D only, #cases > #features)|' ...
               'Linear Local Tangent Space Alignment           (low-D only, >32 GB RAM for high-D data)|' ...
               'Large-Margin Nearest Neighbour                 ()|' ...
               'Linear Discriminant Analysis                   (low-D only)|' ...
               'Neighborhood Component Analysis                (low-D only)'];
    menuact = { 'PCA', ...
                'RobPCA',... 
                'NNMF', ...
                'optNMF', ...
                'NeNMF', ...
                'PLS', ...
                'SparsePCA', ...
                'ProbPCA', ...
                'AutoEncoder', ...
                'FactorAnalysis',...
                'LPP',...
                'LLTSA',...
                'LMNN',...
                'LDA',...
                'NCA' };    
    if isfield(DR,'RedMode'), def = find(strcmp(menuact,DR.RedMode)); else def = 1; end
    if EXPERT
        menustr = [menustr '|' ...
               '* General Discriminant Analysis                (low-D only)|' ...
               '* Multidimensional Scaling                     (low-D only)|' ...
               '* Fast Maximum Variance Unfolding              ()|' ...
               '* Laplacian Eigenmaps                          ()|' ...
               '* Local Linear Embedding                       ()|' ...
               '** Hessian LLE                                 (low-D only, slow, but more accurate than LLE)|' ...
               '** t-Stochastic Nearest Neighbor Embedding     (results highly variable)|' ...
               '** Symmetric SNE                               (results highly variable)|' ...
               '** Sammon''s mapping                            (Sammmon)|', ...
               '** Diffusion Maps                              (slow)|'...
               '** Gaussian Process Latent Variable Modelling  (slow)'];
           menuact = [ menuact, ...
                'GDA',...
                'MDS',...
                'FastMVU',...
                'Laplacian',...
                'LLE',...
                'HessianLLE',...
                't-SNE',...
                'SymSNE',...
                'Sammon',...
                'DiffMaps',...
                'GPLVM'];
    end
    nk_PrintLogo
    mestr = 'Dimensionality reduction techniques'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    if isempty(def), def=1; end
    if EXPERT
        fprintf('\n'); cprintf('red','*  No reconstruction of data in original space supported!');
        fprintf('\n'); cprintf('*red','** No reconstruction of data in original space supported! Only estimate of out-of-sample projection available!');
    end
    act = char(nk_input(mestr,0,'mq', menustr, menuact, def));
    if ~strcmp(act,'BACK'), DR.RedMode = act; end
    if ~exist('PX','var'), PX = []; end

    switch act

        case 'PCA'
            DR.DRsoft = 1; PX = nk_AddParam([], [], [], PX,'reset');
            
        case 'SparsePCA'
            DR.DRsoft = 1;
            if isfield(DR,'SparsePCA')
                lambda = DR.SparsePCA.lambda;
                stop = DR.SparsePCA.stop;
                maxiter = DR.SparsePCA.maxiter;
            else
                lambda = Inf;
                stop = [ .1 .2 .3 .4 .5 ];
                maxiter = 300;
            end
            DR.SparsePCA.lambda     =   nk_input('Lambda parameter (positive ridge term coefficient; set to Inf for high-D problems)',0,'e', lambda);
            DR.SparsePCA.stop       =   nk_input('Stop parameter (if positive => upper bound on the L1-norm of the BETA coefficients, if 0 => regular PCA)',0,'e', stop);
            DR.SparsePCA.maxiter    =   nk_input('max # of iterations)',0,'w1', maxiter);
            PX = nk_AddParam(DR.SparsePCA.lambda, 'SparsePCA-Lambda', 1, PX, 'replace');
            PX = nk_AddParam(DR.SparsePCA.stop, 'SparsePCA-Stop', 1, PX);
            PX = nk_AddParam(DR.SparsePCA.maxiter, 'SparsePCA-Iter', 1, PX);
            
        case {'RobPCA','FactorAnalysis','LDA','GDA','MDS'}
            PX = nk_AddParam([], [], [], PX,'reset');
            
        case 'PLS'
            PX = nk_AddParam([], [], [], PX,'reset');
            uselabel                = nk_input('Use NM label(s) for covariance computation',0,'m','NM label(s)|External label matrix',[1,2],1);
            switch uselabel 
                case 1
                    DR.PLS.behav_mat = nk_MakeDummyVariables(NM.label);
                case 2
                    DR.PLS.behav_mat = nk_input('Provide a side label matrix for covariance matrix computation',0,'e',[],[size(NM.label,1) Inf]);
            end
%             DR.PLS.algostr = char(nk_input('Choose PLS method',0,'m','PLS|Sparse PLS',{'pls','spls'}));
%             switch DR.PLS.algostr
%                 case 
%             
        case 'NNMF'
            
            PX = nk_AddParam([], [], [], PX,'reset');
            NMFmethod               = {'nmf','sparsenmf','convexnmf','orthnmf'};
            NMFmethoddesc           = {'Standard Sparse Non-Negative Matrix Factorization',...
                                      'Sparse Non-Negative Matrix Factorization', ...
                                      'Convex Non-Negative Matrix Factorization', ...
                                      'Orthogonal Non-Negative Matrix T-Factorization'};
            if ~isfield(DR,'NMFmethod'), DR.NMFmethod = NMFmethod{1}; end
            DR.NMFmethod            = NMFmethod{nk_input('Choose NMF method',0,'m',strjoin(NMFmethoddesc,'|'),1:numel(NMFmethoddesc), find(strcmp(NMFmethod,DR.NMFmethod)))};
            DR.options.dis          = 0;
            
        case 'ProbPCA'
            if isfield(DR,'ProbPCA'), iter= DR.ProbPCA.iter; else iter = 200; end
            DR.ProbPCA.iter = nk_input('# of iterations ',0,'e',iter);
            PX = nk_AddParam(DR.ProbPCA.iter, 'ProbPCA-Iter', 1, PX, 'replace');

        case 'LPP'
            if isfield(DR,'LPP'), 
                Sigma = DR.LPP.Sigma; 
                K = DR.LPP.k; 
            else
                Sigma = 1; K = 12;
            end
            DR.LPP.Sigma = nk_input('Bandwith of the Gaussian kernel (Sigma)',0,'e',Sigma);
            DR.LPP.k = nk_input('# of nearest neighbors (k)',0,'e',K);
            PX = nk_AddParam(DR.LPP.Sigma, 'LPP-Sigma', 1, PX, 'replace');
            PX = nk_AddParam(DR.LPP.k, 'LPP-K', 1, PX);
            %if isfield(DR,'LPP') && isfield(DR.LPP,'constuctW'), options = DR.LPP.constructW.options; else options = []; end
    %         [DR.LPP.constructW.options, PX] = nk_ConstructW_config(options, [], navistr);
    %         
    %         fprintf('\n\n *********** Linear Graph Embedding (LGE) Setup ***********\n\n')
    %         DR.LPP.options.regu = nk_input('LGE Setup: Regularization',0,'m', 'a* = argmax (a''X''WXa)/(a''X''DXa+ReguAlpha*I)|PCA',[1,0],2); 
    %         if DR.LPP.options.regu
    %             DR.LPP.options.ReguAlpha = nk_input('LGE Setup -> ReguAlpha',0,'e',0.1);  
    %             PX = nk_AddParam(DR.LPP.options.ReguAlpha, 'LPP-RegulAlpha', 1, PX);
    %             DR.LPP.options.ReguType = 'Ridge';
            %end
            %DR.DRsoft = 1
            
        case 'AutoEncoder'
            if isfield(DR,'AutoEncoder'), Lambda = DR.AutoEncoder.Lambda; else Lambda = 1; end
            DR.AutoEncoder.Lambda = nk_input('Lambda (L2-regularization parameter)',0,'e',Lambda);
            PX = nk_AddParam(DR.AutoEncoder.Lambda, 'AutoEncoder-Lambda', 1, PX, 'replace');

        case 'LLTSA'
            if isfield(DR,'LLTSA'), k=DR.LLTSA.k; else k=12; end
            DR.LLTSA.k = nk_input('# of nearest neighbors (k)',0,'e',k);                 
            PX = nk_AddParam(DR.LLTSA.k, 'LLTSA-K', 1, PX, 'replace');

        case 'Laplacian'
            if isfield(DR,'Laplacian'), k=DR.Laplacian.k; else k=12; end
            DR.Laplacian.k = nk_input('# of nearest neighbors (k)',0,'e',k);                 
            PX = nk_AddParam(DR.Laplacian.k, 'Laplacian-K', 1, PX, 'replace');
            DR.Laplacian.s = nk_input('Sigma',0,'e',1);              
            PX = nk_AddParam(DR.Laplacian.s, 'Laplacian-Sigma', 1, PX);

        case 'LLE'
            if isfield(DR,'LLE'), k = DR.LLE.k; else k=12; end
            DR.LLE.k = nk_input('# of nearest neighbors (k)',0,'e',k);
            PX = nk_AddParam(DR.LLE.k, 'LLE-K', 1, PX, 'replace');

        case 'HessianLLE'
            if isfield(DR,'HessianLLE'), k = DR.HessianLLE.k; else k=12; end
            DR.HessianLLE.k = nk_input('k = ',0,'e',k); 
            PX = nk_AddParam(DR.HessianLLE.k, 'HessianLLE-K', 1, PX, 'replace');

        case 't-SNE'
            if isfield(DR,'tSNE'), 
                dims = DR.tSNE.PCAdims; perplex = DR.tSNE.perplex;
                if isfield(DR.tSNE,'PCAPercMode'), PercMode = DR.tSNE.PCAPercMode; end
            else
                dims = []; perplex = 30; PercMode = []; 
            end
            [DR.tSNE.PCAdims, DR.tSNE.PCAPercMode] = return_dims('Initial PCA', PercMode, dims);
            DR.tSNE.perplex = nk_input('t-SNE perplexity [5-50]',0,'e',perplex);
            PX = nk_AddParam(DR.tSNE.PCAdims, 'tSNE-PCA-dimensions',  1, PX, 'replace');
            PX = nk_AddParam(DR.tSNE.perplex, 'tSNE-perplexity', 1, PX);

        case 'SymSNE'
            if isfield(DR,'SymSNE'), 
                perplex = DR.SymSNE.perplex;
            else
                perplex = 30;
            end
            DR.SymSNE.perplex = nk_input('SymSNE perplexity [5-50]',0,'e',perplex);  
            PX = nk_AddParam(DR.SymSNE.perplex, 'SymSNE-perplexity', 1, PX, 'replace');

        case 'FastMVU'
            if isfield(DR,'FastMVU'), 
                k = DR.FastMVU.k; finetune = DR.FastMVU.finetune; if ~finetune; finetune = 2; end
            else
                k = 12; finetune = 0;
            end
            DR.FastMVU.k = nk_input('# of nearest neighbors (k)',0,'e',k);                    
            PX = nk_AddParam(DR.k, 'FastMVU-K', 1, PX, 'replace');
            DR.finetune = nk_input('Fine-Tuning (only out-of sample estimation supported in this mode)?',0,'yes|no',[1,0], finetune);
        
        case 'NCA'
            DR.NCA.lambda = nk_input('Lambda range',0,'e',0.5);
            PX = nk_AddParam(DR.NCA.lambda, 'NCA-Lambda', 1, PX, 'replace');

        case 'Sammon'
            PX = nk_AddParam([], [], [], PX,'reset');

        case 'DiffMaps'
            if isfield(DR,'DiffMaps'), Sigma = DR.DiffMaps.Sigma; else Sigma = 1; end
            if isfield(DR,'DiffMaps'), T = DR.DiffMaps.t; else T = 1; end
            DR.DiffMaps.Sigma = nk_input('Variance of the Gaussian (Sigma) for the affinity computation',0,'e',Sigma);
            DR.DiffMaps.t = nk_input('Operator applied on the graph (t)',0,'e',T);
            PX = nk_AddParam(DR.DiffMaps.Sigma, 'DiffMaps-Sigma', 1, PX, 'replace');
            PX = nk_AddParam(DR.DiffMaps.t, 'DiffMaps-T', 1, PX);

        case 'GPLVM'
            if isfield(DR,'GPLVM'), Sigma = DR.GPLVM.Sigma; else Sigma = 1; end
            DR.GPLVM.Sigma= nk_input('Filter bandwith (Sigma)',0,'e',Sigma);
            PX = nk_AddParam(DR.GPLVM.Sigma, 'GPLVM-Sigma', 1, PX, 'replace');

        case 'LMNN'
            PX = nk_AddParam([], [], [], PX,'reset');
            %DR.k = nk_input('K neighbors (currently uses only default = 3 NN)',0,'e',3);
            %[PX.Params, PX.Params_desc] = nk_AddParam(DR.k, 'K', PX.Params, PX.Params_desc);
    end
else
    DR.RedMode = 'PCA'; DR.DRsoft = 1; PX = nk_AddParam([], [], [], PX,'reset');
end

end