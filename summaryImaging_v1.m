function summaryImaging_v1
clear
cd('F:\Imaging in GC\RVKC368\180731')
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
clearvars -except data

cd('F:\Imaging in GC\RVKC377\180821')
load('data.mat')
data2 = resp3;
for i = 1:length(data2)
    data2(i).trace_dF = neuron(i).trace_dF;
    data2(i).framT = neuron(i).framT;
    data2(i).Licktrace_dF = neuron(i).Licktrace_dF;
    data2(i).TLick        = neuron(i).TLick;
    data2(i).T = neuron(i).T;
    data2(i).S_Taste_dF = neuron(i).S_Taste_dF;
    data2(i).M_Taste_dF = neuron(i).M_Taste_dF;
    data2(i).CA_Taste_dF = neuron(i).CA_Taste_dF;
    data2(i).Q_Taste_dF = neuron(i).Q_Taste_dF;
    data2(i).W_Taste_dF = neuron(i).W_Taste_dF;
end
Session1 = [data,data2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Session1
cue = [Session1.CueRes];
lick = [Session1.LickRes];
prop.cue = sum(cue)./length(cue);
licknum = find([Session1.CueRes]==0 & [Session1.LickRes] ==1);
prop.lick = length(licknum)/length(Session1);
for i = 1:length(Session1)
    Taste(i).S = Session1(i).Sres;
    Taste(i).M = Session1(i).Mres;
    Taste(i).CA = Session1(i).CAres;
    Taste(i).Q = Session1(i).Qres;
    Taste(i).W = Session1(i).Wres;
end

Taste = cell2mat(squeeze(struct2cell(Taste))');
TasteRes = find(sum(Taste,2)>=1);
TasteRes = setdiff(TasteRes,licknum);
prop.taste = length(TasteRes)/length(Session1);

save('Session1.mat','prop','cue','licknum','Taste','TasteRes','Session1')

%% Session 2
clear
cd('F:\Imaging in GC\RVKC368\180802')
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
clearvars -except data

cd('F:\Imaging in GC\RVKC377\180823')
load('data.mat')
data2 = resp3;
for i = 1:length(data2)
    data2(i).trace_dF = neuron(i).trace_dF;
    data2(i).framT = neuron(i).framT;
    data2(i).Licktrace_dF = neuron(i).Licktrace_dF;
    data2(i).TLick        = neuron(i).TLick;
    data2(i).T = neuron(i).T;
    data2(i).S_Taste_dF = neuron(i).S_Taste_dF;
    data2(i).M_Taste_dF = neuron(i).M_Taste_dF;
    data2(i).CA_Taste_dF = neuron(i).CA_Taste_dF;
    data2(i).Q_Taste_dF = neuron(i).Q_Taste_dF;
    data2(i).W_Taste_dF = neuron(i).W_Taste_dF;
end
Session2 = [data,data2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Session2
cue = [Session2.CueRes];
lick = [Session2.LickRes];
prop.cue = sum(cue)./length(cue);
licknum = find([Session2.CueRes]==0 & [Session2.LickRes] ==1);
prop.lick = length(licknum)/length(Session2);
for i = 1:length(Session2)
    Taste(i).S = Session2(i).Sres;
    Taste(i).M = Session2(i).Mres;
    Taste(i).CA = Session2(i).CAres;
    Taste(i).Q = Session2(i).Qres;
    Taste(i).W = Session2(i).Wres;
end

Taste = cell2mat(squeeze(struct2cell(Taste))');
TasteRes = find(sum(Taste,2)>=1);
TasteRes = setdiff(TasteRes,licknum);
prop.taste = length(TasteRes)/length(Session2);

cd('F:\Imaging in GC\Summary')
save('Session2.mat','prop','cue','licknum','Taste','TasteRes','Session2')

%% Session 3
clear
cd('F:\Imaging in GC\RVKC368\180807')
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
clearvars -except data

cd('F:\Imaging in GC\RVKC377\180828')
load('data.mat')
data2 = resp3;
for i = 1:length(data2)
    data2(i).trace_dF = neuron(i).trace_dF;
    data2(i).framT = neuron(i).framT;
    data2(i).Licktrace_dF = neuron(i).Licktrace_dF;
    data2(i).TLick        = neuron(i).TLick;
    data2(i).T = neuron(i).T;
    data2(i).S_Taste_dF = neuron(i).S_Taste_dF;
    data2(i).M_Taste_dF = neuron(i).M_Taste_dF;
    data2(i).CA_Taste_dF = neuron(i).CA_Taste_dF;
    data2(i).Q_Taste_dF = neuron(i).Q_Taste_dF;
    data2(i).W_Taste_dF = neuron(i).W_Taste_dF;
end
Session3 = [data,data2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Session3
cue = [Session3.CueRes];
lick = [Session3.LickRes];
prop.cue = sum(cue)./length(cue);
licknum = find([Session3.CueRes]==0 & [Session3.LickRes] ==1);
prop.lick = length(licknum)/length(Session3);
for i = 1:length(Session3)
    Taste(i).S = Session3(i).Sres;
    Taste(i).M = Session3(i).Mres;
    Taste(i).CA = Session3(i).CAres;
    Taste(i).Q = Session3(i).Qres;
    Taste(i).W = Session3(i).Wres;
end

Taste = cell2mat(squeeze(struct2cell(Taste))');
TasteRes = find(sum(Taste,2)>=1);
TasteRes = setdiff(TasteRes,licknum);
prop.taste = length(TasteRes)/length(Session3);
cd('F:\Imaging in GC\Summary')
save('Session3.mat','prop','cue','licknum','Taste','TasteRes','Session3')