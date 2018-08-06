function resp = tasteResponse(neuron)
t = 0.05;
rw = 3;
for j = 1:length(neuron)
    % j = 15;
    idx = find(neuron(j).T>-1 & neuron(j).T <0);
    T_idx1 = find(neuron(j).T>0 & neuron(j).T <rw);
%     T_idx2 = find(neuron(j).T>1 & neuron(j).T <2);
%     T_idx3 = find(neuron(j).T>2 & neuron(j).T <3);
    S_baseline    = mean(neuron(j).S_Taste_dF(:,idx),2);
    S_Taste_1     = mean(neuron(j).S_Taste_dF(:,T_idx1),2);
%     S_Taste_2   = mean(neuron(j).S_Taste_dF(:,T_idx2),2);
%     S_Taste_3    = mean(neuron(j).S_Taste_dF(:,T_idx3),2);
    
    [p(1),h(1)] = signrank(S_Taste_1-S_baseline);
    if mean(S_Taste_1)< mean(S_baseline) || p(1)>t || mean(S_Taste_1)<0;
        h(1) = 0;
    end
%     [p(2),h(2)] = ranksum(S_baseline,S_Taste_2);
%     if mean(S_Taste_2)< mean(S_baseline);
%         h(2) = 0;
%     end
%     [p(3),h(3)] = ranksum(S_baseline,S_Taste_3);
%     if mean(S_Taste_3)< mean(S_baseline);
%         h(3) = 0;
%     end
    M_baseline    = mean(neuron(j).M_Taste_dF(:,idx),2); % 2nd taste
    M_Taste_1     = mean(neuron(j).M_Taste_dF(:,T_idx1),2);
    [p(2),h(2)] = signrank(M_Taste_1-M_baseline);
    if mean(M_Taste_1)< mean(M_baseline)|| p(2)>t || mean(M_Taste_1)<0 ;
        h(2) = 0;
    end
    
    CA_baseline    = mean(neuron(j).CA_Taste_dF(:,idx),2); % 3rd taste
    CA_Taste_1     = mean(neuron(j).CA_Taste_dF(:,T_idx1),2);
    [p(3),h(3)] = signrank(CA_Taste_1-CA_baseline);
    if mean(CA_Taste_1)< mean(CA_baseline)|| p(3)>t;
        h(3) = 0;
    end
    
    Q_baseline    = mean(neuron(j).Q_Taste_dF(:,idx),2); % 4th taste
    Q_Taste_1     = mean(neuron(j).Q_Taste_dF(:,T_idx1),2);
    [p(4),h(4)] = signrank(Q_Taste_1-Q_baseline);
    if mean(Q_Taste_1)< mean(Q_baseline)|| p(4)>t;
        h(4) = 0;
    end   
    
    W_baseline    = mean(neuron(j).W_Taste_dF(:,idx),2); % 4th taste
    W_Taste_1     = mean(neuron(j).W_Taste_dF(:,T_idx1),2);
    [p(5),h(5)] = signrank(W_Taste_1-W_baseline);
    if mean(W_Taste_1)< mean(W_baseline)|| p(5)>t || mean(W_Taste_1)<0;
        h(5) = 0;
    end    
    resp(j).Sres = h(1);
    resp(j).Mres = h(2);
    resp(j).CAres= h(3);
    resp(j).Qres = h(4);
    resp(j).Wres = h(5);
end
%% cue response
t = 0.05;
rw = 2;
for j = 1:length(neuron)
    % j = 15;
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).framT>0 & trial(1).framT <rw);
%     C_idx2 = find(trial(1).framT>1 & trial(1).framT <2);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).traceSmooth_dF(j,C_idx1),2);
%         C_2(i)        = mean(trial(i).traceSmooth_dF(j,C_idx1),2);
    end
    [c_p(1),c_h(1)] = signrank(C_1-baseline);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 ;
        c_h(1) = 0;
    end
%     [c_p(2),c_h(2)] = ranksum(baseline,C_2);
%     if mean(C_2)< mean(baseline) || c_p(2)>t || mean(C_2)<0 ;;
%         c_h(2) = 0;
%     end
    resp(j).CueRes = c_h(1);
%     neuron(j).CueRes.h = sum(h);
%     neuron(j).TasteResponse = sum(h);
end
%% statistical test here the baseline is the 1 s before the cue; test lick response
t = 0.05;
rw = 1;
for j = 1:length(neuron)
    % j = 15;
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).TLick>0 & trial(1).TLick <rw);
%     C_idx2 = find(trial(1).framT>1 & trial(1).framT <2);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).LickTrace(j,C_idx1),2);
%         C_2(i)        = mean(trial(i).traceSmooth_dF(j,C_idx1),2);
    end
    [c_p(1),c_h(1)] = signrank(C_1-baseline);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 ;
        c_h(1) = 0;
    end
%     [c_p(2),c_h(2)] = ranksum(baseline,C_2);
%     if mean(C_2)< mean(baseline) || c_p(2)>t || mean(C_2)<0 ;;
%         c_h(2) = 0;
%     end
    resp(j).LickRes = c_h(1);
%     neuron(j).CueRes.h = sum(h);
%     neuron(j).TasteResponse = sum(h);
end
