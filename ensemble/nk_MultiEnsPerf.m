function [hEPerf, hE] = nk_MultiEnsPerf(E, sE, L, C)
global MULTI RAND

switch MULTI.method
    case 1 % Simple One-Vs-One / One-vs-All
        [hE, hEPerf] = nk_MultiDecideMaxWins(E, L, C, MULTI.decisiontype, 1);
    case 2 % Error-Correcting Output Codes
        if ~isfield(RAND,'Decompose'),
            decompose = 1;
        else
            decompose = RAND.Decompose;
        end
        [hE, hEPerf] = nk_MultiDecideErrorCorrOutCodes(sE, L, C, decompose, MULTI.decoding, 1);
    case 3 % Hierarchical One-Vs-One
        [hE, hEPerf] = nk_MultiDecideHierarchOneVsOne(E, L, C, MULTI.decisiontype);
    case 4 % Directed Acyclic Graph
        [hE, hEPerf] = nk_MultiDecideDAG(E, L, C, MULTI.decisiontype, 1);
%     case 5
%         [hE, hEPerf] = nk_MultiDecideMaxWinsNegTiesOpt(E, L, C, MULTI.decisiontype, 1);
end

return