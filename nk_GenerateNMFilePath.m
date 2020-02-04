function pth = nk_GenerateNMFilePath(RootDir,AnalDesc,DatType, OptDesc, StrOut, DatID, CV2Perm , CV2Fold, CV1Perm, CV1Fold)

CV2Desc = sprintf('_oCV%g.%g', CV2Perm, CV2Fold);
if exist('CV1Perm','var') && exist('CV1Fold','var') && ~isempty(CV1Perm) && ~isempty(CV1Fold) 
    CV2Desc = sprintf('_iCV%g.%g', CV1Perm, CV1Fold);
end
pth = fullfile(RootDir, sprintf('%s_%s%s%s%s_ID%s.mat', AnalDesc, DatType, OptDesc, CV2Desc, StrOut, DatID)); 

