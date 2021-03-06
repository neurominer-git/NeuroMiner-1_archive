function load_selModality(handles)

for i=1:numel(handles.visdata)
    popuplist{i} = sprintf('Modality %g: %s', i, handles.NM.datadescriptor{handles.visdata{i}.params.varind}.desc);
end
set(handles.selModality, 'String', popuplist);

popuplist=[];
if isfield(handles.visdata{1},'MEAN'),                      popuplist{1} = 'Feature weights [Overall Mean (StErr)]';                        end
if isfield(handles.visdata{1},'MEAN_CV2'),                  popuplist{end+1} = 'Feature weights [Grand Mean (StErr)]';                      end
if isfield(handles.visdata{1},'CVRatio'),                   popuplist{end+1} = 'CV-ratio of feature weights [Overall Mean]';                end
if isfield(handles.visdata{1},'CVRatio_CV2'),               popuplist{end+1} = 'CV-ratio of feature weights [Grand Mean]';                  end
if isfield(handles.visdata{1},'FeatProb'),                  popuplist{end+1} = 'Feature selection probability [Overall Mean]';              end
if isfield(handles.visdata{1},'Prob_CV2'),                  popuplist{end+1} = 'Probability of feature reliability (95%-CI) [Grand Mean]';  end
if isfield(handles.visdata{1},'Spearman_CV2'),              popuplist{end+1} = 'Spearman correlation [Grand Mean]';                         end
if isfield(handles.visdata{1},'Pearson_CV2'),               popuplist{end+1} = 'Pearson correlation [Grand Mean]';                          end
if isfield(handles.visdata{1},'Spearman_CV2_p_uncorr'),     popuplist{end+1} = 'Spearman correlation -log10(P value) [Grand Mean]';         end
if isfield(handles.visdata{1},'Pearson_CV2_p_uncorr'),      popuplist{end+1} = 'Pearson correlation -log10(P value) [Grand Mean]';          end
if isfield(handles.visdata{1},'Spearman_CV2_p_fdr'),        popuplist{end+1} = 'Spearman correlation -log10(P value, FDR) [Grand Mean]';    end
if isfield(handles.visdata{1},'Pearson_CV2_p_fdr'),         popuplist{end+1} = 'Pearson correlation -log10(P value, FDR) [Grand Mean]';     end
if isfield(handles.visdata{1},'PermProb_CV2'),              popuplist{end+1} = 'Permutation-based -log10(P value) [Grand Mean]';            end
if isfield(handles.visdata{1},'PermProb_CV2_FDR_PVAL'),     popuplist{end+1} = 'Permutation-based -log10(P value, FDR) [Grand Mean]';       end
if isfield(handles.visdata{1},'PermZ_CV2'),                 popuplist{end+1} = 'Permutation-based Z Score [Grand Mean]';                    end
if isfield(handles.visdata{1},'PermModel_Eval_Global'),     popuplist{end+1} = 'Model P value histogram';                                   end
set(handles.selVisMeas, 'String', popuplist); 