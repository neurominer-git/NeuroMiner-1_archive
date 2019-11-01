function [CURACT, act ] = nk_PartialCorrelations_config(NM, CURACT, varind, parentstr, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
% Defaults
COVAR_DEF       = 1;
COVDIR_DEF      = 1;
INTERCEPT_DEF   = 2;
BETAEXT_DEF     = [];
SUBGROUP_DEF    = [];

if ~defaultsfl
    
    % Get information from CURACT if available
    if isfield(CURACT,'COVAR'),     COVAR_DEF       = CURACT.COVAR; end
    if isfield(CURACT,'INTERCEPT'), INTERCEPT_DEF   = CURACT.INTERCEPT; end
    if isfield(CURACT,'COVDIR'),    COVDIR_DEF      = CURACT.COVDIR; end
    if isfield(CURACT,'SUBGROUP'),  SUBGROUP_DEF    = CURACT.SUBGROUP; end
    if isfield(CURACT,'BETAEXT'),   BETAEXT_DEF     = CURACT.BETAEXT; end
    
    COVAR_STR = sprintf('%s, ', NM.covnames{COVAR_DEF}); COVAR_STR(end-1:end)=[];
    if INTERCEPT_DEF == 2,      INTERCEPT_STR = 'yes'; else     INTERCEPT_STR = 'no'; end; menuact = [1 2];
    if COVDIR_DEF == 1,         COVDIR_STR = 'attenuate'; else  COVDIR_STR = 'increase'; end; menuact = [menuact 3];
    SUBGROUP_MNU1=[]; SUBGROUP_MNU2=[];
    if ~isempty(BETAEXT_DEF),   
        BETAEXT_STR = 'yes';
        if isfinite(BETAEXT_DEF), 
            BETAEXT_MAT = sprintf('%g x %g matrix defined', size(BETAEXT_DEF,1), size(BETAEXT_DEF,2)); 
        else
            BETAEXT_MAT = 'undefined'; 
        end
        BETAEXT_MNU = sprintf('|Define beta coeficients [ %s ]',BETAEXT_MAT);
        menuact = [ menuact 4 5 ];
        
    else
        BETAEXT_STR = 'no'; 
        BETAEXT_MNU = [];
        menuact = [ menuact 4 ];
        if ~isempty(SUBGROUP_DEF), 
            SUBGROUP_STR = 'yes'; 
            if isfinite(SUBGROUP_DEF), 
                SUBGROUP_MAT = sprintf('vector with %g case(s) defined', sum(SUBGROUP_DEF)); 
            else
                SUBGROUP_MAT = 'undefined';
            end
            SUBGROUP_MNU2 = sprintf('|Provide index to training cases for computing betas [ %s ]', SUBGROUP_MAT );
            menuact = [menuact 6 7];
        else
            SUBGROUP_STR = 'no'; SUBGROUP_MNU2 = [];
            menuact = [menuact 6];
        end
        SUBGROUP_MNU1 = sprintf('|Define subgroup of training cases for computing betas [ %s ]',  SUBGROUP_STR);
    end
    
    menustr = ['Select covariates [ ' COVAR_STR ' ]', ...
               '|Include intercept [ ' INTERCEPT_STR ' ]', ...
               '|Attenuate or Increase covariate effects [ ' COVDIR_STR ' ]' ...
               '|Use externally-computed beta coeficients [ ' BETAEXT_STR ' ]' ...
               BETAEXT_MNU ...
               SUBGROUP_MNU1 ...
               SUBGROUP_MNU2];
    
    nk_PrintLogo
    mestr = 'Residualization setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            COVAR_DEF = nk_SelectCovariateIndex(NM, COVAR_DEF, 1);
   
        case 2
            if INTERCEPT_DEF == 2, INTERCEPT_DEF = 1; elseif INTERCEPT_DEF == 1, INTERCEPT_DEF = 2; end
        
        case 3
            if COVDIR_DEF == 1, COVDIR_DEF = 2; elseif COVDIR_DEF == 2, COVDIR_DEF = 1; end
        
        case 4
            if ~isfield(CURACT,'BETAEXT'), 
                CURACT.BETAEXT = NaN; 
            elseif isfield(CURACT,'BETAEXT'), 
                CURACT = rmfield(CURACT,'BETAEXT'); 
            end
        case 5
            if INTERCEPT_DEF
                defc = numel(COVAR_DEF) + 1;
            else
                defc = numel(COVAR_DEF) ;
            end
            CURACT.BETAEXT = nk_input('Define precompute beta matrix',0,'e',[],[defc, size(NM.Y{varind},2)]);

        case 6
              if ~isfield(CURACT,'SUBGROUP'), 
                CURACT.SUBGROUP = NaN; 
            elseif isfield(CURACT,'SUBGROUP'), 
                CURACT = rmfield(CURACT,'SUBGROUP'); 
              end            
        case 7
            CURACT.SUBGROUP = logical(nk_input('Define index vector (ones = used / zeros = unused) for beta computation',0,'e',[],[numel(NM.label),1]));
    end
end

CURACT.COVAR        = COVAR_DEF;
CURACT.INTERCEPT    = INTERCEPT_DEF;
CURACT.COVDIR       = COVDIR_DEF;
