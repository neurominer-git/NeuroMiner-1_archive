function [ Tr, CV, Ts, Ocv, TrainedParam ] = nk_ReturnAtOptPos(oTr, oCV, oTs, oOcv, Pnt, z)

Ocv = [];
if ~isempty(Pnt.data_ind)
    Ix = Pnt.data_ind(z);
else
    Ix = 1;
end
Tr = oTr{Ix}; CV = oCV{Ix}; Ts = oTs{Ix};   
if ~isempty(oOcv), Ocv = oOcv{Ix}; end

if ~isempty(Pnt.nA)
    TrainedParam = cell(1,Pnt.nA);
    for a = 1:Pnt.nA
        if isstruct(Pnt.TrainedParam{a})
            TrainedParam{a} = Pnt.TrainedParam{a};
        else
            TrainedParam{a} = Pnt.TrainedParam{a}(Pnt.train_ind(z,a));
        end
    end
else
    TrainedParam = Pnt.TrainedParam;
end
