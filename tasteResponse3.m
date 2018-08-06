function resp = tasteResponse2(neuron,trial)
t = 0.05;
rw = 3;
for j = 1:length(neuron)
    % j = 15;
    idx = find(neuron(j).T>-1 & neuron(j).T <0);
    T_idx1 = find(neuron(j).T>0 & neuron(j).T <1);
    T_idx2 = find(neuron(j).T>1 & neuron(j).T <2);
    T_idx3 = find(neuron(j).T>2 & neuron(j).T <3);
    S_baseline    = mean(neuron(j).S_Taste_dF(:,idx),2);
    S_Taste_1     = mean(neuron(j).S_Taste_dF(:,T_idx1),2);
    S_Taste_2     = mean(neuron(j).S_Taste_dF(:,T_idx2),2);
    S_Taste_3     = mean(neuron(j).S_Taste_dF(:,T_idx3),2);
    S_Taste       = [S_Taste_1, S_Taste_2, S_Taste_3];
    for i = 1:size(S_Taste,2)
        [p(i),h(i)] = signrank(S_Taste(:,i)-S_baseline);
        if p(i)>0.05/3 || mean(S_Taste(:,i))<0 || mean(S_Taste(:,i))< mean(S_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        s_resp = 1;
    else
        s_resp = 0;
    end
   
    M_baseline    = mean(neuron(j).M_Taste_dF(:,idx),2); % 2nd taste
    M_Taste_1     = mean(neuron(j).M_Taste_dF(:,T_idx1),2);
    M_Taste_2     = mean(neuron(j).M_Taste_dF(:,T_idx2),2);
    M_Taste_3     = mean(neuron(j).M_Taste_dF(:,T_idx3),2);
    M_Taste       = [M_Taste_1, M_Taste_2, M_Taste_3];
    for i = 1:size(M_Taste,2)
        [p(i),h(i)] = signrank(M_Taste(:,i)-M_baseline);
        if p(i)>0.05/3 || mean(M_Taste(:,i))<0 || mean(M_Taste(:,i))< mean(M_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        m_resp = 1;
    else
        m_resp = 0;
    end
    
    CA_baseline    = mean(neuron(j).CA_Taste_dF(:,idx),2); % 3rd taste
    CA_Taste_1     = mean(neuron(j).CA_Taste_dF(:,T_idx1),2);
    CA_Taste_2     = mean(neuron(j).CA_Taste_dF(:,T_idx2),2);
    CA_Taste_3     = mean(neuron(j).CA_Taste_dF(:,T_idx3),2);
    CA_Taste       = [CA_Taste_1, CA_Taste_2, CA_Taste_3];
    for i = 1:size(CA_Taste,2)
        [p(i),h(i)] = signrank(CA_Taste(:,i)-CA_baseline);
        if p(i)>0.05/3 || mean(CA_Taste(:,i))<0 || mean(CA_Taste(:,i))< mean(CA_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        ca_resp = 1;
    else
        ca_resp = 0;
    end
    
    Q_baseline    = mean(neuron(j).Q_Taste_dF(:,idx),2); % 4th taste
    Q_Taste_1     = mean(neuron(j).Q_Taste_dF(:,T_idx1),2);
    Q_Taste_2     = mean(neuron(j).Q_Taste_dF(:,T_idx2),2);
    Q_Taste_3     = mean(neuron(j).Q_Taste_dF(:,T_idx3),2);
    Q_Taste       = [Q_Taste_1, Q_Taste_2, Q_Taste_3];
    for i = 1:size(Q_Taste,2)
        [p(i),h(i)] = signrank(Q_Taste(:,i)-Q_baseline);
        if p(i)>0.05/3 || mean(Q_Taste(:,i))<0 || mean(Q_Taste(:,i))< mean(Q_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        q_resp = 1;
    else
        q_resp = 0;
    end   
    
    W_baseline    = mean(neuron(j).W_Taste_dF(:,idx),2); % 4th taste
    W_Taste_1     = mean(neuron(j).W_Taste_dF(:,T_idx1),2);
    W_Taste_2     = mean(neuron(j).W_Taste_dF(:,T_idx2),2);
    W_Taste_3     = mean(neuron(j).W_Taste_dF(:,T_idx3),2);
    W_Taste       = [W_Taste_1, W_Taste_2, W_Taste_3];
    for i = 1:size(W_Taste,2)
        [p(i),h(i)] = signrank(W_Taste(:,i)-W_baseline);
        if p(i)>0.05/3 || mean(W_Taste(:,i))<0 || mean(W_Taste(:,i))< mean(W_baseline) 
            h(i) = 0;
        end
    end
    if sum(h)>0
        w_resp = 1;
    else
        w_resp = 0;
    end  
    resp(j).Sres = s_resp;
    resp(j).Mres = m_resp;
    resp(j).CAres = ca_resp;
    resp(j).Qres = q_resp;
    resp(j).Wres = w_resp;
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
        [p(i),h(i)] = signrank(C(:,i)-baseline);
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
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).TLick>0 & trial(1).TLick <rw);
    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).LickTrace(j,C_idx1),2);
    end
    [c_p(1),c_h(1)] = signrank(C_1-baseline);
    if mean(C_1)< mean(baseline) || c_p(1)>t || mean(C_1)<0 ;
        c_h(1) = 0;
    end

    resp(j).LickRes = c_h(1);

end

