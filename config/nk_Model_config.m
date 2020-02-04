function [SVM, act] = nk_Model_config(SVM, TrainParam, parentstr, varind)

global NM EXPERT

d = nk_GetParamDescription2(NM,TrainParam,'GridParam');
if (isfield(NM,'TrainParam') && isfield(NM.TrainParam,'RAND') && NM.TrainParam.RAND.InnerFold == -1 )
    SVM.GridParam = 1; optimstr = ''; menusel = [];
else
    optimstr = d.GridParam; optimstr = ['Define ML model performance criterion [ ' optimstr ' ]|']; menusel = 3;
end

if isfield(SVM,'prog')
    d = nk_GetParamDescription2(NM,TrainParam,'SVMprog',d);
    d = nk_GetParamDescription2(NM,TrainParam,'classifier',d);
    d = nk_GetParamDescription2(NM,TrainParam,'kernel',d);
    if isempty(strfind(d.prog,'RVM'))
        progstr = [d.prog ' => ' d.classifier];
    else
        progstr = d.prog;
    end
    progstr = [ 'Configure learning algorithm [ ' progstr ' ]|'];
    kernelstr = ['Define kernel type [ ' d.kernel ' ]|']; 
    
    % Some learners do not operate in kernel space - skip kernel menu item
    % in these cases
    if any(strcmp(SVM.prog,{'IMRELF','kNNMEX','LIBLIN','GLMFIT','MEXELM', 'DECTRE','RNDFOR','GLMNET','GRDBST','SEQOPT'}))
        mnuact = [optimstr progstr ];
        mnusel = [menusel 1];
    elseif strcmp(SVM.prog,'matLRN')
        cprintf('blue','\nDefine possible kernel configuration required for a matLearn algorithm in the Learning algorithm parameters setup')
        mnuact = [optimstr progstr ];
        mnusel = [menusel 1];
    else
        mnuact = [optimstr progstr kernelstr] ;
        mnusel = [menusel 1 2];
    end
else
    
    if isempty(SVM)
        cprintf('red','No model configuration found for modality %g', varind)
        copyflag = nk_input('Do you want to use another modality''s configuration as template for current model configuration',0,'yes|no',[1,0],1);
        if copyflag
            varind_copy = nk_input('Specify modality index',0,'e');
            if varind_copy <= numel(NM.TrainParam.SVM) && ~isempty(NM.TrainParam.SVM{varind_copy})
                SVM = NM.TrainParam.SVM{varind_copy};
            end
        end
    end
    
    progstr ='NA';
    optimstr = 'NA';
    mnuact = [optimstr progstr];
    mnusel = [3 1];
end

if EXPERT && ~strcmp(SVM.prog, 'SEQOPT')
    if isfield(SVM,'Post') && isfield(SVM.Post,'Detrend') && SVM.Post.Detrend
        detrendstr = 'enabled';
    else
        detrendstr = 'Disabled / NA';
    end
    switch NM.modeflag
        case 'regression'
            mnuact = [ mnuact 'Detrend predicted targets [ ' detrendstr ']|' ]; 
            mnusel = [ mnusel 4 ];

        case 'classification'
            mnuact = [ mnuact 'Optimize decision threshold based on ROC [ ' detrendstr ' ]|' ]; 
            mnusel = [ mnusel 5 ];
    end
end

if isfield(SVM,'kernel'), kerntype = SVM.kernel.kernstr; else kerntype = []; end

switch NM.modeflag
    case 'classification'
        if EXPERT 
            adasyndef = {'yes','no'};
            if ~isfield(SVM,'ADASYN'),
                SVM.ADASYN.flag = 2;
                adasynstr = adasyndef{2}; 
            else
                adasynstr = adasyndef{SVM.ADASYN.flag};
            end    
            mnuact = [ mnuact 'Use ADASYN to adjust for unbalanced class settings [ ' adasynstr ' ]|' ];
            mnusel = [ mnusel 6 ];
            if SVM.ADASYN.flag == 1
                if ~isfield(SVM.ADASYN,'beta'), SVM.ADASYN.beta = 1; end
                betadef = SVM.ADASYN.beta; 
                mnuact = [ mnuact sprintf('Define beta value (defines how much balancing will be applied, 0<->1) [ k=%g ]|',betadef)];
                mnusel = [ mnusel 7 ];
                if ~isfield(SVM.ADASYN,'kDensity'), SVM.ADASYN.kDensity = 5; end
                kdensdef = SVM.ADASYN.kDensity; 
                mnuact = [ mnuact sprintf('Define k value of density algorithm (kNN looking at both classes) [ k=%g ]|',kdensdef)];
                mnusel = [ mnusel 8 ];
                if ~isfield(SVM.ADASYN,'kSMOTE'), SVM.ADASYN.kSMOTE = 5; end
                ksmotedef = SVM.ADASYN.kSMOTE; 
                mnuact = [ mnuact sprintf('Define k value of SMOTE algorithm (kNN looking only at the minority class) [ k=%g ]|',ksmotedef)];
                mnusel = [ mnusel 9 ];
                if ~isfield(SVM.ADASYN,'normalized'), SVM.ADASYN.normalized = false; end
                if SVM.ADASYN.normalized, normstr = 'yes'; else, normstr = 'no'; end
                mnuact = [ mnuact sprintf('Data is already normalized for kNN search [ %s ]|',normstr) ];
                mnusel = [ mnusel 10 ];
            end
        end
    case 'regression'
        SVM.ADASYN.flag = 2;
end

switch NM.modeflag
    
    case 'classification'
        if EXPERT
            sftmenu = ['SVM --------------> LIBSVM|' ...
                       'SVM / LR ---------> LIBLINEAR|' ...
                       'ROBUBST SVM ------> Outlier-insensitive LIBSVM|' ...
                       'NEURAL -----------> Extreme Learning Machine (you have to scale data to -1/1 and prune low-var or NaN vectors) |' ...
                       'RVM / BAYES ------> Psorakis & Damoulas''s implementation with Multiple Kernel Learning|' ...
                       'RVM / BAYES ------> Relevance vector classification (MC_RVM)|' ...
                       'RVM / BAYES ------> Bayesian Sparse Logistic Regression (for classification in high-D spaces)|' ...
                       'LOCAL LEARNING ---> k-Nearest Neighbors (kNN)|' ...
                       'LOCAL LEARNING ---> IMRelief|' ...
                       'UNIVARIATE -------> GLM Logistic Regression|' ...
                       'DECISION TREE ----> MATLAB Statistics Toolbox implementation of the Decision Tree algorithm|' ...
                       'RANDOM FOREST ----> Abhishek Jaiantilal''s fast Random Forest algorithm for MATLAB|' ...
                       'matLearn ---------> Mark Schmidt''s matLearn library|'...
                       'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ...
                       'GRADIENT BOOST----> Carlos Becker''s gradient boosting algorithm'];
           sftval   = [ 'LIBSVM'; ...
                        'LIBLIN'; ...
                        'ROBSVM'; ...
                        'MEXELM'; ...
                        'MKLRVM'; ...
                        'MVTRVR'; ...
                        'BLOREG'; ...
                        'kNNMEX'; ...
                        'IMRELF'; ...
                        'GLMFIT'; ...
                        'DECTRE'; ...
                        'RNDFOR'; ...
                        'matLRN'; ...
                        'GLMNET'; ...
                        'GRDBST'];
            if TrainParam.STACKING.flag==1
                sftmenu = [sftmenu ...
                        '|SEQOPT -----------> Predictive sequence optimization algorithm (Stacking only)'];
                sftval = [sftval; ...
                        'SEQOPT'];
            end
        else
            sftmenu = ['SVM --------------> LIBSVM|' ...
                       'SVM/LR -----------> LIBLINEAR|' ...
                       'LOCAL LEARNING ---> k-Nearest Neighbors (kNN)|' ...
                       'RVM / BAYES ------> Relevance vector classification (MC_RVM)|' ...
                       'UNIVARIATE -------> GLM Logistic Regression|' ...
                       'DECISION TREE ----> MATLAB Statistics Toolbox implementation of the Decision Tree algorithm|' ...
                       'RANDOM FOREST ----> Abhishek Jaiantilal''s fast Random Forest algorithm for MATLAB|'...
                       'matLearn ---------> Mark Schmidt''s matLearn library|'...
                       'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ...
                       'GRADIENT BOOST----> Carlos Becker''s gradient boosting algorithm'];
           
            sftval   = ['LIBSVM'; ...
                        'LIBLIN'; ...
                        'kNNMEX'; ...
                        'MVTRVR'; ...
                        'GLMFIT'; ...
                        'DECTRE'; ...
                        'RNDFOR'; ...
                        'matLRN'; ...
                        'GLMNET'; ...
                        'GRDBST'];
        end
              
    case 'regression'
        if EXPERT
            sftmenu = [ 'SVM --------------> LIBSVM|' ...
                        'SVM --------------> LIBLINEAR|' ...
                        'RVM / BAYES ------> Mike Tipping''s algorithm|' ...
                        'RVM / BAYES ------> Multinomial relevance vector regression (MVRVR) (for multiple target label prediction)|' ...
                        'UNIVARIATE -------> GLM Linear Regression|' ...
                        'RANDOM FOREST ----> Abhishek Jaiantilal''s fast Random Forest algorithm for MATLAB|' ...
                        'matLearn ---------> Select a machine learning strategy from Mark Schmidt''s matLearn toolbox|' ...
                        'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ... 
                        'GRADIENT BOOST----> Carlos Becker''s gradient boosting algorithm'];
            sftval = ['LIBSVM'; ...
                      'LIBLIN'; ...
                      'MikRVM'; ...
                      'MVTRVR'; ...
                      'GLMFIT'; ...
                      'RNDFOR'; ...
                      'matLRN'; ...
                      'GLMNET'; ...
                      'GRDBST'];
        else
            sftmenu = [ 'SVM --------------> LIBSVM|' ...
                        'SVM --------------> LIBLINEAR|' ...
                        'RVM / BAYES ------> Mike Tipping''s algorithm|' ...
                        'UNIVARIATE -------> GLM Linear Regression|' ...
                        'RANDOM FOREST ----> Abhishek Jaiantilal''s fast Random Forest algorithm for MATLAB|' ...
                        'matLearn ---------> Select a machine learning strategy from Mark Schmidt''s matLearn toolbox|' ...
                        'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ... 
                        'GRADIENT BOOST----> Carlos Becker''s gradient boosting algorithm'];
            sftval = ['LIBSVM'; ...
                      'LIBLIN'; ...
                      'MikRVM'; ...
                      'GLMFIT'; ...
                      'RNDFOR'; ...
                      'matLRN'; ...
                      'GLMNET'; ...
                      'GRDBST'];
        end
end

nk_PrintLogo
mestr = 'ML algorithm parameters'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr);
act = nk_input(mestr, 0,'mq', mnuact, mnusel,1);

switch act
 
    case 1
        nk_PrintLogo
        mestr = ['Define model algorithm for ' NM.modeflag]; cprintf('*blue','\nYou are here: %s >>>',navistr);
        selProg = nk_input(mestr,0,'mq', sftmenu, 1:size(sftval,1), 1);
        if selProg, 
            
            SVM.prog = sftval(selProg,:);
            
            switch SVM.prog

                case {'LIBSVM','LIBLIN'}
                   
                    SVM = nk_SVM_config(NM, SVM, SVM.prog, kerntype, navistr);
                    SVM = nk_Kernel_config(SVM);
                
                case 'MikRVM'
                    SVM.RVMflag = true;
                    SVM.RVM.UserOpt = SB2_UserOptions;
                    SVM.RVM.ParamSet = SB2_ParameterSettings;
                    switch NM.modeflag
                        case 'regression'
                           sftval = char(nk_input('Select likelihood model',0,'m', ...
                                    ['Gaussian (for real valued functions)|' ...
                                     'Poisson (for positive integer value function)'],[1 2],1));
                    end
                    if sftval
                        sftlst = {'GAUSS','POISS'};
                        SVM.RVM.LikelihoodModel = sftlst{sftval};
                    end
                    SVM = nk_Kernel_config(SVM);
                    
                case 'MKLRVM'
                    SVM = nk_MKLRVM_config(SVM);
                    SVM = nk_Kernel_config(SVM);

                case 'IMRELF'
                    SVM.imrelief = nk_IMRelief_config(SVM);
                    SVM.kernel.kerndesc = 'other';
                    SVM.kernel.kerndef = 1;
                    SVM.kernel.kernstr = ' -t 2';

                case 'GLMFIT'
                    SVM.RVMflag = true; % Probability flag
                    SVM = nk_Kernel_config(SVM);

                case 'kNNMEX'
                    SVM.RVMflag = true;
                    
                case 'KPCSVM'
                    SVM = nk_KPCA_LINSVM(SVM);

                case 'BLOREG'
                    SVM = nk_Kernel_config(SVM);

                case 'LSTSVM'
                    SVM = nk_LSTSVM_config(SVM);
                    
                case 'MVTRVR'
                    SVM.MVTRVR.iter = nk_input('# of iterations (Less # give sparser, but less precise models)',0,'e',20);
                    SVM = nk_Kernel_config(SVM);
                
                case 'FAMRVR'
                    SVM.FAMRVR.iter = nk_input('# of iterations of the EM algorithm',0,'e',1000);
                    SVM.FAMRVR.tolerance = nk_input('Tolerance',0,'e',.1);
                    SVM = nk_Kernel_config(SVM);
                    
                case {'MEXELM','DECTRE', 'RNDFOR'}
                    %SVM.MEXELM.nHidden = nk_input('# of hidden neurons',0,'e',100);
                    switch SVM.prog
                        case 'DECTRE'
                            if ~isfield(SVM,SVM.prog), SVM.(SVM.prog) = []; end
                        otherwise
                            SVM = nk_Kernel_config(SVM);
                    end
                    
                case 'matLRN'
                    if ~isfield(SVM,'matLRN'), SVM.matLRN = []; end
                    switch NM.modeflag
                        case 'classification'
                            matLRN = nk_matLearn_config(SVM.matLRN,'binaryclass',1);
                           
                        otherwise
                            matLRN = nk_matLearn_config(SVM.matLRN,'regression',1);
                    end
                    if isempty(SVM.matLRN) || ~strcmp(matLRN.algo{1},SVM.matLRN.algo{1})
                        NM.TrainParam.GRD.matLearn = nk_matLearn_config(matLRN,matLRN.learner.framework,3);
                    end
                    SVM.matLRN = matLRN;
                    
                case {'GLMNET','GRDBST'}
                    if ~isfield(SVM,SVM.prog), SVM.(SVM.prog) = []; end
                    if ~isfield(NM.TrainParam.GRD,SVM.prog), NM.TrainParam.GRD.(SVM.prog) = nk_GLMNET_config(SVM.prog, [],1); end
                    switch SVM.prog
                        case 'GLMNET'
                            switch NM.modeflag
                                case 'classification'
                                    SVM.GLMNET.family = 'binomial';
                                otherwise
                                    SVM.GLMNET.family = 'gaussian';
                            end
                    end
                    
                case 'ROBSVM'
                    if ~isfield(SVM,SVM.prog), SVM.(SVM.prog) = []; end
                    if ~isfield(NM.TrainParam.GRD,SVM.prog), NM.TrainParam.GRD.(SVM.prog) = nk_ROBSVM_config(SVM.prog, [], NM.modeflag, 1); end
                    
                case 'SEQOPT'
                    act = 1; if ~isfield(SVM,'SEQOPT'), [~, SVM.SEQOPT ] = nk_SEQOPT_config([],1); end
                    while act 
                      [act, SVM.SEQOPT ] = nk_SEQOPT_config( SVM.SEQOPT);
%                     SVM.(SVM.prog).C = nk_input('Enter candidate sequences to be compared against each other (matrix padded with NaNs)',0,'e');
%                     SVM.(SVM.prog).Lims = nk_input('Enter minimum lower and maximum upper percentiles for ambiguous case forwarding',0,'e');
%                     SVM.(SVM.prog).nCutOff = 5;                    
                    end
            end
        end

    case 2  
        
        SVM = nk_Kernel_config(SVM);
        
    case 3

        SVM.GridParam = nk_EvalFunc_config(NM, SVM, navistr);
        
    case 4
        
        SVM.Post.Detrend = nk_input('Enable post-hoc linear prediction detrending (using C1 test data)',0,'yes|no',[1,0],1);
        
    case 5
        
        SVM.Post.Detrend = nk_input('Enable post-hoc optimization of decision threshold (using ROC analysis of C1 test data)',0,'yes|no',[1,0],1);
        
    case 6
        
        if SVM.ADASYN.flag ==1, SVM.ADASYN.flag = 2; else, SVM.ADASYN.flag = 1; end
        
    case 7
        
        SVM.ADASYN.beta = nk_input('Define beta for degree of balancing',0,'e',SVM.ADASYN.beta);
        
    case 8
        
        SVM.ADASYN.kDensity = nk_input('Define k for Density estimation in ADASYN',0,'i',SVM.ADASYN.kDensity);
        
    case 9
        
        SVM.ADASYN.kSMOTE = nk_input('Define k for SMOTE in ADASYN',0,'i',SVM.ADASYN.kSMOTE);
        
    case 10
        
        SVM.ADASYN.normalized = ~SVM.ADASYN.normalized;
        
    
end

