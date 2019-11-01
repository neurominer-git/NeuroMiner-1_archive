function P = nk_CalibrateProbabilities(P)

% indN = P<0;
% indP = P>0;
% P(indN) = P(indN)+0.5;
% P(indP) = P(indP)-0.5;
P = P-0.5;

end