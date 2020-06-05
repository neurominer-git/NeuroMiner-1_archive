% =========================================================================
% FORMAT NM = nk_DefineOOCVData2_config(NM, O, parentstr)
% =========================================================================
% This function coordinates the input of independent test data into NM
%
% Inputs:
% -------
% NM        : The NM workspace
% O         : The OOCV workspace (containt
% parentstr : Name of the calling function
%
% Outputs:
% --------
% NM (see above)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 09/2017
function NM = nk_DefineOOCVData2_config(NM, O, parentstr)
    
nM = numel(NM.Y);

if (~exist('O','var') || isempty(O)) || isnumeric(O)
    [NM, Y, oocvind, fldnam, dattype] = nk_SelectOOCVdata(NM, O, 1);
else

    fldnam      = O.fldnam;
    oocvind     = O.ind;
    Y           = NM.(fldnam){oocvind};
    Y.desc      = O.desc;
    Y.date      = O.date;
    Y.label_known = false;
    
    switch O.fldnam
        case 'OOCV'
            dattype     = 'independent test data';
        case 'C'
            dattype     = 'calibration data';
    end
end
if isempty(oocvind), return; end
NM.(fldnam){oocvind}.desc = Y.desc;
NM.(fldnam){oocvind}.date = Y.date;
if ~isfield(NM.(fldnam){oocvind},'n_subjects_all'), NM.(fldnam){oocvind}.n_subjects_all = Inf; end
if ~isfield(NM.(fldnam){oocvind},'labels_known') && isfield(Y,'labels_known') || (NM.(fldnam){oocvind}.labels_known ~=  Y.labels_known)
    NM.(fldnam){oocvind}.labels_known = Y.labels_known; 
end

na_str = '?';

if isfield(NM,'covars'), covflag = true; else, covflag = false; end

% Loop through modalities
for i=1:nM

    nk_PrintLogo
    fprintf('\n\n'); mestr = sprintf('Input %s for Modality %g', dattype, i);  navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','\nYou are here: %s >>>',navistr); 
    
    % Retrieve input settings from the discovery data
    IO = NM.datadescriptor{i}.input_settings;
    if isfield(IO,'selCases'), IO = rmfield(IO,'selCases'); end
    
    % Remove setting for non-labeled subjects
    IO.nangroup=false;
    IO.nan_subjects=0;
    if isfield(IO,'Pnan')
        IO = rmfield(IO,'Pnan');
        IO = rmfield(IO,'Vnan');
    end

    % Activate independent test data input
    IO.oocvflag = true;
    IO.labels_known = Y.labels_known;
    IO.badcoords = NM.badcoords{i}; 
    IO.brainmask = NM.brainmask{i}; 
    IO.Ydims = size(NM.Y{i},2);
    
    if IO.labels_known
        IO.n_subjects = IO.n_subjects/0;
        IO.n_subjects_all = Inf;
    else
        IO.n_subjects = Inf;
        IO.n_subjects_all = Inf;
        IO.n_samples = 1;
    end
    
    if i>1 && isfield(NM.(fldnam){oocvind},'cases')
        IO.ID = NM.(fldnam){oocvind}.cases;
    else
        IO = rmfield(IO,'ID');
        if isfield(IO,'survanal_time'), IO = rmfield(IO,'survanal_time'); end
    end
    
    if strcmp(IO.datasource,'matrix')
        IO.matrix_edit = na_str;
        IO.sheets = na_str;
        IO.sheet = na_str;
        IO.sheets = na_str;
        IO.M_edit = na_str;
        IO.featnames_cv = NM.featnames{i};
    else
        if strcmp(IO.datasource,'spm')
            IO.datasource = 'nifti'; 
            IO.groupmode = 1;
            IO = rmfield(IO,'design');
            IO = SetFileFilter(IO,IO.groupmode,IO.datasource);
        end
        IO.globvar_edit = na_str;
        if isfield(IO,'g') && ~isempty(IO.g), IO = rmfield(IO,'g'); end
        if IO.labels_known
            IO.P = repmat({[]},1,IO.n_samples);
            IO.V = IO.P;
        else
            IO.P = []; 
            IO.V = [];
        end
        IO.PP = [];
        IO = rmfield(IO,'Vinfo');
        IO = rmfield(IO,'Vvox');
        IO = rmfield(IO,'F');
        IO = rmfield(IO,'files');
        if isfield(IO,'L') && ~isempty(IO.L),
            IO.label_edit = na_str;
            IO = rmfield(IO,'L'); 
        end
    end
    t_act = Inf; t_mess = [];while ~strcmp(t_act,'BACK'), [ IO, t_act, t_mess ] = DataIO( NM.(fldnam)(oocvind) , mestr, IO, t_mess, i);  end
    if IO.completed
        NM.(fldnam){oocvind} = TransferModality2NM( NM.(fldnam){oocvind}, IO, i ); 
        NM.(fldnam){oocvind}.n_subjects_all = size(NM.(fldnam){oocvind}.cases,1);
    else
        NM.(fldnam){oocvind} = Y;
        break;
    end
end

% Don't forget the covariates if they are present in the discovery data
if covflag && isfield(NM.(fldnam){oocvind},'Y') && NM.(fldnam){oocvind}.n_subjects_all>0
   NM.(fldnam){oocvind}.covars = nk_DefineCovars_config(NM.(fldnam){oocvind}.n_subjects_all, NM.covars); 
end