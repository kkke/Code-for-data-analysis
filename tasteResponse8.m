function resp = tasteResponse8(neuron,trial)
t = 0.05;
rw = 3.5;
taste = {'S','M','CA','Q','W'};
for j = 1:length(neuron)
    % j = 15;
    idx = find(neuron(j).T>-1 & neuron(j).T <0);
    T_idx1 = find(neuron(j).T>0 & neuron(j).T <rw);
    for i = 1:length(taste)
        baseline    = mean(neuron(j).([taste{i},'_Taste_dF'])(:,idx),2); % 4th taste
        TasteR     = mean(neuron(j).([taste{i},'_Taste_dF'])(:,T_idx1),2);
        [p,h] = ranksum(TasteR,baseline);
        if p>=t
            h = 0;
        elseif p<t && mean(TasteR)>mean(baseline)
            h = 1;
        elseif p<t && mean(TasteR)<mean(baseline)
            h = -1;
        end
        resp(j).([taste{i},'res']) = h;
    end
    %     resp(j).Sres = s_resp;
%     resp(j).Mres = m_resp;
%     resp(j).CAres = ca_resp;
%     resp(j).Qres = q_resp;
%     resp(j).Wres = w_resp;
end
%% cue response with signed rank
clear p h
t = 0.01;
rw = 2;
for j = 1:length(neuron)
    idx = find(neuron(1).framT>-1 & neuron(1).framT <0); % use [-2,-1] as the baseline
    C_idx1 = find(neuron(1).framT>0 & neuron(1).framT <2);
    baseline   = mean(neuron(j).trace_dF(:,idx),2);
    C_1        = mean(neuron(j).trace_dF(:,C_idx1),2);

    C = C_1;
    [p,h] = ranksum(C,baseline);
    if p >=t
        h = 0;
    elseif p<t && mean(C)< mean(baseline)
        h = -1;
    elseif p<t && mean(C)> mean(baseline)
        h = 1;
    end
    resp(j).CueRes = h;
%     resp(j).cueP = p;
%     neuron(j).CueRes.h = sum(h);
%     neuron(j).TasteResponse = sum(h);
end
%% lick response with signed ranksum test
clear baseline C_1
t = 0.01;
rw = 1;
for j = 1:length(neuron)
    idx = find(trial(1).framT>-1 & trial(1).framT <0); % [-2,-1] as the baseline
    C_idx1 = find(trial(1).TLick>0 & trial(1).TLick <rw);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).LickTrace(j,C_idx1),2);
    end
    [c_p,c_h] = ranksum(C_1,baseline);
    if  c_p>=t
        c_h = 0;
    elseif c_p<t && mean(C_1)< mean(baseline)
        c_h = -1;
    elseif c_p<t && mean(C_1)> mean(baseline)
        c_h = 1;
    end

    resp(j).LickRes = c_h;

end

