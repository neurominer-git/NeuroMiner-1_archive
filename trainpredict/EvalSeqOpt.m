function [CritGain, ExamFreq, PercThreshU, PercThreshL, AbsThreshU, AbsThreshL] = EvalSeqOpt(Models)

[ix, jx, nclass] = size(Models);
CritGain = zeros(ix,jx,nclass);
ExamFreq = cell(ix,jx,nclass);
PercThreshU = ExamFreq;
PercThreshL = ExamFreq;
AbsThreshU = ExamFreq;
AbsThreshL = ExamFreq;
for h=1:nclass
    for k = 1:ix
        for l= 1:jx
            CritGain(k,l,h) = Models{k,l}{h}.(Models{k,l}{h}.Crit)(end) - Models{k,l}{h}.(Models{k,l}{h}.Crit)(1); 
            ExamFreq{k,l,h} = Models{k,l}{h}.examsfreq; 
            PercThreshU{k,l,h} = Models{k,l}{h}.optuvec; 
            PercThreshL{k,l,h} = Models{k,l}{h}.optlvec;
            AbsThreshU{k,l,h} = Models{k,l}{h}.optuthr;
            AbsThreshL{k,l,h} = Models{k,l}{h}.optlthr;
        end
    end
end