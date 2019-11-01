function [act, SEQOPT ] = nk_SEQOPT_config(SEQOPT, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end

if ~defaultsfl 
    
    MODEARR = {'Only population to be propagated','Entire population'}; MODESTR = MODEARR{SEQOPT.Mode};
    MnuStr = sprintf('Define training population for optimization [ %s ]', MODESTR);                                                    MnuAct = 1;
    
    LIMSSTR = sprintf('Lower: %g%%, Upper: %g%%', SEQOPT.Lims(1), SEQOPT.Lims(2));
    MnuStr = sprintf('%s|Define minimum lower and maximum upper percentiles for ambiguous case propagation [ %s ]',MnuStr, LIMSSTR);    MnuAct = [MnuAct 2];
    
    ANCHORARR = {'Decision boundary','Median'}; ANCHORSTR = ANCHORARR{SEQOPT.AnchorType};
    MnuStr = sprintf('%s|Define of propagation algorithm''s anchor [ %s ]',MnuStr, ANCHORSTR);                                          MnuAct = [MnuAct 3];
    
    nk_PrintLogo
    act = nk_input('Select action',0,'mq', ...
                        MnuStr, ...
                        MnuAct, 1);
    switch act
        case 1
            SEQOPT.Mode = nk_input('Define training population flag for optimization',0,'mq','Only population to be propagated|Entire population',[1,2],SEQOPT.Mode);
        case 2
            SEQOPT.Lims = nk_input('Enter minimum lower and maximum upper percentiles for ambiguous case propagation',0,'e',SEQOPT.Lims);
        case 3
            SEQOPT.AnchorType = nk_input('Anchor stepping to model''s decision boundary or decision score median',0,'mq','Decision boundary|Median',[1,2], SEQOPT.AnchorType);
    end
else
    act =0; SEQOPT = struct('Mode', 1, 'Lims', [0 100], 'AnchorType', 1); 
end


