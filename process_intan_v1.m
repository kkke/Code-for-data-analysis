%% process data recorded from Intan board
function [data,trial] = process_intan_v1(filename);

data = read_Intan(filename);
% data = read_intan_batch;
%% extract the event data
[Frame,~] = Timing_onset_offset(data.event(7,:), data.ts, 0.5,10,1);
[Tone,~]  = Timing_onset_offset(data.event(6,:), data.ts, 0.5,10,1);
[S,~]     = Timing_onset_offset(data.event(1,:), data.ts, 0.5,10,1);
% close all
%% re-organize the data in a trial format
unit  = spike2eventRasteandPSTH_NP(Frame,Tone,100,-3000,12000);
unit2 = spike2eventRasteandPSTH_NP(S,Tone,100,-3000,12000);
for i = 1:length(Tone)
    trial(i).tone = Tone(i);
    trial(i).Frame = unit.spikeraster(i).times;
    if isempty(unit2.spikeraster(i).times)
        trial(i).S = NaN;
    else
        trial(i).S = unit2.spikeraster(i).times
    end
end
%% load the data from scope
% a             = textread('RVKC340_051218.txt');
% [licks,~]     = Timing_onset_offset(a(:,4), a(:,1), -0.22,30,1);
% [paw,~]       = Timing_onset_offset(a(:,4), a(:,1), -0.02,30,1);
% [Tone_sco,~]  = Timing_onset_offset(a(:,2), a(:,1), 2,30,1);
% licks         = setdiff(licks,paw);
% unit2 = spike2eventRasteandPSTH_NP(licks,Tone_sco,100,-2000,12000);
%%
% figure
% for i = 1:length(Tone)
%     if ~isempty(unit.spikeraster(i).times)
%         scatter(unit.spikeraster(i).times,i*ones(size(unit.spikeraster(i).times)),'r+')
%         hold on
%     end
%     if ~isempty(unit2.spikeraster(i).times)
%         scatter(unit2.spikeraster(i).times,i*ones(size(unit2.spikeraster(i).times)),6,'kv','filled')
%         hold on
%     end
% end
% xlim([-2,12])
