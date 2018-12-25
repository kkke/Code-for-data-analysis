%% load all date files
file = {'RVKC368','RVKC377','RVKC402','RVKC403','RVKC404','RVKC405'};
date(1).date = {'180731','180802','180807','180809','180814'};
date(2).date = {'180821','180823','180828'};
date(3).date = {'181022','181024','181029'};
date(4).date = {'181022','181024','181029','181030','181031'};
date(5).date = {'181112','181114','181119','181120'};
date(6).date = {'181112','181114','181119','181120'};
for i = 1:length(file)
    data(i).mice = file{i};
    data(i).date = date(i).date;
end
%% automatically segment all imaging files
for i = 1:length(file)
    for j = 1:length(data(i).date)
        fprintf(['Processing ', data(i).mice, ' ', data(i).date{j}, '\n'])
        segmentation_Suite2P(data(i).mice,data(i).date{j})
    end
end
%% automatically process all segmented data
for i = 1:length(file)
    for j = 1:length(data(i).date)
        fprintf(['Processing ', data(i).mice, ' ', data(i).date{j}, '\n'])
        imaging_analysis_GC_v4(data(i).mice,data(i).date{j})
    end
end
%%


