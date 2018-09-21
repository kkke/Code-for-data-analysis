function resp = tasteResponse4(session)
t = 0.05;
rw = 3;
neuron = session;
taste = {'S_Taste_dF','M_Taste_dF','CA_Taste_dF','Q_Taste_dF','W_Taste_dF'};
respID  = {'Sres','Mres','Cres','Qres','Wres'};
for j = 1:length(neuron)
    for i = 1:length(taste)
        idx = find(neuron(j).T>-1 & neuron(j).T <0);
        T_idx1 = find(neuron(j).T>0 & neuron(j).T <1);
        T_idx2 = find(neuron(j).T>1 & neuron(j).T <2);
        T_idx3 = find(neuron(j).T>2 & neuron(j).T <3);
        baseline    = mean(neuron(j).(taste{i})(:,idx),2);
        Taste_1     = mean(neuron(j).(taste{i})(:,T_idx1),2);
        Taste_2     = mean(neuron(j).(taste{i})(:,T_idx2),2);
        Taste_3     = mean(neuron(j).(taste{i})(:,T_idx3),2);
        Taste       = [Taste_1, Taste_2, Taste_3];
        for k = 1:size(Taste,2)
            [p(k),h(k)] = signrank(Taste(:,k),baseline);
            if p(k)>0.05/3 || mean(Taste(:,k))<0 || mean(Taste(:,k))< mean(baseline)
                h(k) = 0;
            end
        end
        if sum(h)>0
            s_resp = 1;
        else
            s_resp = 0;
        end
        resp(j).(respID{i}) = s_resp;
        clear h
        
    end
end
%% cue response with signed rank
clear p h
t = 0.05;
rw = 2;
for j = 1:length(neuron)
    idx = find(neuron(1).framT>-1 & neuron(1).framT <0);
    C_idx1 = find(neuron(1).framT>0 & neuron(1).framT <1);
    C_idx2 = find(neuron(1).framT>1 & neuron(1).framT <2);
    baseline   = mean(neuron(j).trace_dF(:,idx),2);
    C_1        = mean(neuron(j).trace_dF(:,C_idx1),2);
    C_2        = mean(neuron(j).trace_dF(:,C_idx2),2);

    C = [C_1,C_2];
    for i = 1:size(C,2)
        [p(i),h(i)] = signrank(C(:,i),baseline);
        if p(i)>0.05/2 || mean(C(:,i))<0 || mean(C(:,i))< mean(baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        cue_resp = 1;
    else
        cue_resp = 0;
    end 

    resp(j).CueRes = cue_resp(1);
%     neuron(j).CueRes.h = sum(h);
%     neuron(j).TasteResponse = sum(h);
end
%% lick response with signed ranksum test
clear baseline C_1
t = 0.05;
rw = 1;
for j = 1:length(neuron)
    idx = find(neuron(j).framT>-1 & neuron(j).framT <0);
    C_idx1 = find(neuron(j).TLick>0 & neuron(j).TLick <rw);
    baseline      = mean(neuron(j).trace_dF(:,idx),2);
    C_1           = mean(neuron(j).Licktrace_dF(:,C_idx1),2);
    [c_p(1),c_h(1)] = signrank(C_1,baseline);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 ;
        c_h(1) = 0;
    end
    resp(j).LickRes = c_h(1);

end

