function nk_Utilities
global SPMAVAIL CV NM

menustr = ['Compile C++ files for NM|'...
           'Set root paths of neuroimaging tools (SPM/Freesurfer)|' ...
           'Create PreprocData path master|' ...
           'Create CVdatamat path master|' ...
           'Create CVresults path master|' ...
           'Create VISdatamat path master|' ...
           'Create OptPreprocParam path master|' ...
           'Create OptModel path master'];
menuact = 1:8;

if SPMAVAIL
    menustr = [menustr '|Freesurfer: Downsampling, registration to fsaverage & matrix generation (Linux only)']; menuact = [ menuact 9 ];
end
nk_PrintLogo
act = nk_input('Choose utility function',0,'mq', menustr, menuact, 1);
                 
switch act
    case 0
        return
    case 1
        NMdir = fileparts(which('nm'));
        cfilesdir = [NMdir filesep 'cfiles'];
        nk_make(cfilesdir)
    case 2
        neurominerpath = fileparts(which('neurominer.m'));
        imaging_init_path = fullfile(neurominerpath,'imaging_init.mat');
        nk_ImagingInit(neurominerpath, imaging_init_path, 1);
    case 3
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenPreprocMaster2(NM.id, CV(1), [], inp.rootdir, [], 1, inp.varind, inp.varstr, inp.concatfl);
    case 4
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenCVdataMaster2(NM.id, CV(1), [], inp.rootdir, [], 1, inp.varind, inp.varstr, inp.concatfl);
    case 5
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenCVresultsMaster(NM.id,[],inp.rootdir,[],1); 
    case 6
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenDataMaster(NM.id, 'VISdatamat', CV(1),[],inp.rootdir,[],1);
    case 7
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenDataMaster(NM.id, 'OptPreprocParam', CV(1),[],inp.rootdir,[],1);
    case 8
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenDataMaster(NM.id, 'OptModel', CV(1),[],inp.rootdir,[],1);
    case 9
        [P, Y] = nk_FSreslice;
        assignin('base','P',P);
        assignin('base','Y',Y);
end

nk_Utilities

end

