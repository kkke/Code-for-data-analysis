function summaryImaging_v2
clear
% Load each Imaging data from Session 1
file = {'RVKC368\180731','RVKC377\180821','RVKC402\181022','RVKC403\181022'...
    ,'RVKC404\181112','RVKC405\181112'};
data2 = [];
for j = 1:length(file)
    cd(['F:\Imaging in GC\ImagingData\',file{j}])
    load('data.mat')
    data = resp3;
    for i = 1:length(data)
        data(i).trace_dF = neuron(i).trace_dF;
        data(i).framT = neuron(i).framT;
        data(i).Licktrace_dF = neuron(i).Licktrace_dF;
        data(i).TLick        = neuron(i).TLick;
        data(i).T = neuron(i).T;
        data(i).S_Taste_dF = neuron(i).S_Taste_dF;
        data(i).M_Taste_dF = neuron(i).M_Taste_dF;
        data(i).CA_Taste_dF = neuron(i).CA_Taste_dF;
        data(i).Q_Taste_dF = neuron(i).Q_Taste_dF;
        data(i).W_Taste_dF = neuron(i).W_Taste_dF;
    end
    clearvars -except data data2 file
    data2 = [data2,data];
end
