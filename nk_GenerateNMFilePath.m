function pth = nk_GenerateNMFilePath(RootDir,AnalDesc,DatType, OptDesc, StrOut, DatID, CV2Perm ,CV2Fold)

CVDesc = sprintf('_oCV%g.%g', CV2Perm, CV2Fold);
pth = fullfile(RootDir, sprintf('%s_%s%s%s%s_ID%s.mat', AnalDesc, DatType, OptDesc, CVDesc, StrOut, DatID)); 

end