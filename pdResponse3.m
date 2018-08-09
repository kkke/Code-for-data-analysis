function resp = pdResponse3(neuron,trial)
%% cue response with bootstrap
clear p h
t = 0.01;
rw = 1;
for j = 1:length(neuron)
    idx = find(neuron(1).framT>-1 & neuron(1).framT <0);
    C_idx1 = find(neuron(1).framT>0 & neuron(1).framT <1);
    baseline   = mean(neuron(j).trace_dF(:,idx),2);
    C_1        = mean(neuron(j).trace_dF(:,C_idx1),2);
    pdata = [baseline;C_1];
    mean_dif = mean(C_1)-mean(baseline);
    iteration =10000;
    for k = 1:iteration
        y = datasample(pdata,length(pdata));
        pseudo_baseline = y(1:length(baseline));
        pseudo_Cue     = y(length(baseline)+1:end);
        pseudo_mean_dif(k) = mean(pseudo_Cue)-mean(pseudo_baseline);
    end
    p(j) = (length(find(abs(pseudo_mean_dif)>=abs(mean_dif)))+1)./(iteration+1);
    if p(j)>t || mean(C_1)<0 || mean(C_1)< mean(baseline)
        h = 0;
    else
        h = 1;
    end
    resp(j).CueRes = h;
    %     neuron(j).CueRes.h = sum(h);
    %     neuron(j).TasteResponse = sum(h);
end
%% lick response with signed bootstrap
clear baseline C_1
t = 0.01;
rw = 2;
for j = 1:length(neuron)
    idx = find(trial(1).framT>-1 & trial(1).framT <0);
    C_idx1 = find(trial(1).TLick>0 & trial(1).TLick <1);
    C_idx2 = find(trial(1).TLick>1 & trial(1).TLick <2);

    for i = 1:length(trial)
        baseline(i)   = mean(trial(i).traceSmooth_dF(j,idx),2);
        C_1(i)        = mean(trial(i).LickTrace(j,C_idx1),2);
        C_2(i)        = mean(trial(i).LickTrace(j,C_idx2),2);
    end
    C = [C_1;C_2];
    for i = 1:size(C,1)
        pdata = [baseline,C(i,:)];
        mean_dif = mean(C(i,:))-mean(baseline);
        iteration =10000;
        for k = 1:iteration
            y = datasample(pdata,length(pdata));
            pseudo_baseline = y(1:length(baseline));
            pseudo_Cue     = y(length(baseline)+1:end);
            pseudo_mean_dif(k) = mean(pseudo_Cue)-mean(pseudo_baseline);
        end
        c_p(i) = (length(find(abs(pseudo_mean_dif)>=abs(mean_dif)))+1)./(iteration+1);
%         if c_p(i)>t || mean(C(i,:))<0.05 || mean(C(i,:))< mean(baseline)
        if c_p(i)>t || mean(C(i,:))< mean(baseline)|| mean(C(i,:))<0
            c_h(i) = 0;
        else
            c_h(i) = 1;
       end
    end
    if sum(c_h)>0
        lick_resp = 1;
    else
        lick_resp = 0;
    end
    resp(j).LickRes = lick_resp;

end

