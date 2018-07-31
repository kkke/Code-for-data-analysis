%% behaviral analysis for 2p rig
%% Step 1: load data from Intan
% add directory
addpath('C:\Users\ke-roberto\Documents\MATLAB\Code_Ke\')
%% load data
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
%% extract the events
thr = 2;
[data.centSp.dig,~] = Timing_onset_offset(dataRaw.analog(1,:), dataRaw.ts, thr,30,0); % get the central licks
[Taste_5,Taste_5_2]     = Timing_onset_offset(dataRaw.event(5,:), dataRaw.ts, 0.5,30,0);
[Taste_4,Taste_4_2]     = Timing_onset_offset(dataRaw.event(4,:), dataRaw.ts, 0.5,30,0);
[Taste_3,Taste_3_2]     = Timing_onset_offset(dataRaw.event(3,:), dataRaw.ts, 0.5,30,0);
[Taste_2,Taste_2_2]     = Timing_onset_offset(dataRaw.event(2,:), dataRaw.ts, 0.5,30,0);
[Taste_1,Taste_1_2]     = Timing_onset_offset(dataRaw.event(1,:), dataRaw.ts, 0.5,30,0);
Taste_1 = sort([Taste_1,Taste_1_2]);
Taste_2 = sort([Taste_2,Taste_2_2]);
Taste_3 = sort([Taste_3,Taste_3_2]);
Taste_4 = sort([Taste_4,Taste_4_2]);
Taste_5 = sort([Taste_5,Taste_5_2]);

%% reorganize the data
if isempty(Taste_1)
    Taste_1 =[];
else
    Taste_1(2,:) = 0;
end

if isempty(Taste_1)
    Taste_2 =[];
else
    Taste_2(2,:) = 1;
end

if isempty(Taste_1)
    Taste_3 =[];
else
    Taste_3(2,:) = 2;
end

if isempty(Taste_1)
    Taste_4 =[];
else
    Taste_4(2,:) = 3;
end

if isempty(Taste_1)
    Taste_5 =[];
else
    Taste_5(2,:) = 4;
end
Taste= [Taste_1, Taste_2, Taste_3, Taste_4, Taste_5];
[timestamps1,Idx] = sort(Taste(1,:));
for i = 1:length(Idx)
    data1(i) = Taste(2,Idx(i));
end
%%
a=unique(data1);              % get the id of stimuli
message1=['You train the animal with ',num2str(length(a)),' tastant, Please specify the tastant for each line.'];
message2=[num2str((a(:)+1)')];
uiwait(msgbox({message1, message2}));
x=input('Specify the tastant\n');
% if strcmp(stimuli,'Suc')
%     Suc                        = timestamps1(data1==0);
%     data.Sucrose.dig           = Suc(1:2:length(Suc));
%     [data.Sucrose.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.Sucrose.dig, 100, -5000, 5000);
%     for i=1:length(data.Sucrose.psth_raster.spikeraster)
%         if ~ isempty(data.Sucrose.psth_raster.spikeraster(i).times)
%             for j=1:length(data.Sucrose.psth_raster.spikeraster(i).times)
%                 plot([data.Sucrose.psth_raster.spikeraster(i).times(j),data.Sucrose.psth_raster.spikeraster(i).times(j)],[i,i+1],'k')
%                 hold on
%             end
%         else
%         end
%     end
% end
% figure;
switch length(a)
    case 1
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
        [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        for i=1:length(data.tastant_1.psth_raster.spikeraster)
            h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
            hold on
        end
        legend(h1,x{1})
    case 2
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
        [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        tastant_2                                                =timestamps1(data1==a(2));
        data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
        [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        data.tastant_2.id=x{2};
        for i=1:length(data.tastant_1.psth_raster.spikeraster)
            h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
            hold on
        end
        for j=1:length(data.tastant_2.psth_raster.spikeraster)
            i=i+1;
            h2=scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
            hold on
        end
        legend([h1,h2],{x{1},x{2}})
        
    case 3
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
        [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        tastant_2                                                =timestamps1(data1==a(2));
        data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
        [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
        tastant_3                                                =timestamps1(data1==a(3));
        data.tastant_3.dig                                       =tastant_3(1:2:length(tastant_3 ));
        [data.tastant_3.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_3.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        data.tastant_2.id=x{2};
        data.tastant_3.id=x{3};
        for i=1:length(data.tastant_1.psth_raster.spikeraster)
            h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
            hold on
        end
        for j=1:length(data.tastant_2.psth_raster.spikeraster)
            i=i+1;
            h2=scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
            hold on;
        end
        for j=1:length(data.tastant_3.psth_raster.spikeraster)
            i=i+1;
            h3=scatter(data.tastant_3.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_3.psth_raster.spikeraster(j).times)),6,'bv','filled');
            hold on
        end
        legend([h1,h2,h3],{x{1},x{2},x{3}})
    case 4
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
        [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        tastant_2                                                =timestamps1(data1==a(2));
        data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
        [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
        tastant_3                                                =timestamps1(data1==a(3));
        data.tastant_3.dig                                       =tastant_3(1:2:length(tastant_3 ));
        [data.tastant_3.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_3.dig, 100, -5000, 5000);
        tastant_4                                                =timestamps1(data1==a(4));
        data.tastant_4.dig                                       =tastant_4(1:2:length(tastant_4));
        [data.tastant_4.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_4.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        data.tastant_2.id=x{2};
        data.tastant_3.id=x{3};
        data.tastant_4.id=x{4};
        for i=1:length(data.tastant_1.psth_raster.spikeraster)
            if ~isempty(data.tastant_1.psth_raster.spikeraster(i).times)
                h1=scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_2.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_2.psth_raster.spikeraster(j).times)
                h2=scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_3.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_3.psth_raster.spikeraster(j).times)
                h3 = scatter(data.tastant_3.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_3.psth_raster.spikeraster(j).times)),6,'bv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_4.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_4.psth_raster.spikeraster(j).times)
                h4= scatter(data.tastant_4.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_4.psth_raster.spikeraster(j).times)),6,'rv','filled');
                hold on
            end
        end
        legend([h1,h2,h3,h4],{x{1},x{2},x{3},x{4}})
    case 5
        tastant_1                                                =timestamps1(data1==a(1));
        data.tastant_1.dig                                       =tastant_1(1:2:length(tastant_1));
        [data.tastant_1.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_1.dig, 100, -5000, 5000);
        tastant_2                                                =timestamps1(data1==a(2));
        data.tastant_2.dig                                       =tastant_2(1:2:length(tastant_2 ));
        [data.tastant_2.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_2.dig, 100, -5000, 5000);
        tastant_3                                                =timestamps1(data1==a(3));
        data.tastant_3.dig                                       =tastant_3(1:2:length(tastant_3 ));
        [data.tastant_3.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_3.dig, 100, -5000, 5000);
        tastant_4                                                =timestamps1(data1==a(4));
        data.tastant_4.dig                                       =tastant_4(1:2:length(tastant_4));
        [data.tastant_4.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_4.dig, 100, -5000, 5000);
        tastant_5                                                =timestamps1(data1==a(5));
        data.tastant_5.dig                                       =tastant_5(1:2:length(tastant_5));
        [data.tastant_5.psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.tastant_5.dig, 100, -5000, 5000);
        data.tastant_1.id=x{1};
        data.tastant_2.id=x{2};
        data.tastant_3.id=x{3};
        data.tastant_4.id=x{4};
        data.tastant_5.id=x{5};
        for i=1:length(data.tastant_1.psth_raster.spikeraster)
            if ~isempty(data.tastant_1.psth_raster.spikeraster(i).times)
                h1 =scatter(data.tastant_1.psth_raster.spikeraster(i).times, i*ones(size(data.tastant_1.psth_raster.spikeraster(i).times)),6,'cv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_2.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_2.psth_raster.spikeraster(j).times)
                h2 =scatter(data.tastant_2.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_2.psth_raster.spikeraster(j).times)),6,'mv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_3.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_3.psth_raster.spikeraster(j).times)
                h3 = scatter(data.tastant_3.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_3.psth_raster.spikeraster(j).times)),6,'bv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_4.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_4.psth_raster.spikeraster(j).times)
                h4 = scatter(data.tastant_4.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_4.psth_raster.spikeraster(j).times)),6,'rv','filled');
                hold on
            end
        end
        for j=1:length(data.tastant_5.psth_raster.spikeraster)
            i=i+1;
            if ~isempty(data.tastant_5.psth_raster.spikeraster(j).times)
                h5 = scatter(data.tastant_5.psth_raster.spikeraster(j).times, i*ones(size(data.tastant_5.psth_raster.spikeraster(j).times)),6,'kv','filled');
                hold on
            end
        end
        legend([h1,h2,h3,h4,h5],{x{1},x{2},x{3},x{4},x{5}})
    otherwise
        warning('You have put more than 5 tastant. You need to modify the code to process it')
end
set(gca,'TickDir','out')
ylabel('Trial No.','FontSize',12,'FontWeight','bold')
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
save('data.mat','data')
