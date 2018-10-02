function [cen_duration cen_boutC] = lickbout(cen_L)
ili=diff(cen_L);
b=find(ili>0.5);     % inter lick interval is bigger than 500 ms, which means a new licking bouts
for i=1:length(b)
    if i==1;
        bouts{1}=1:b(1);
    else
        bouts{i}=b(i-1)+1:b(i);
    end
end
bouts{i+1}=b(i)+1:length(cen_L); % all the bouts including less than 3 licks
lickSpont=[];
% find out the bouts which has less than 3 licks
for i=1:length(bouts)
    if length(bouts{i})<3
       randlick=bouts{i};
       bouts{i}=[];
       lickSpont=[lickSpont randlick];
    end  
end
bouts=bouts(~cellfun('isempty',bouts)); % all the real bouts 
for i = 1:length(bouts)
    lickbout{i} = cen_L(bouts{i});
    cen_duration(i) = lickbout{i}(end) - lickbout{i}(1);
    cen_boutC(i)       = length(lickbout{i});
end
figure;
histogram(cen_duration,20)