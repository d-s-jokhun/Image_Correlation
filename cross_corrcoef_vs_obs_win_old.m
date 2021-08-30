clear all

duration=2 %in min
[~,~,raw_pos_xls]=xlsread('rec_4FPS_005_rawtelo.xls');

FPS=(size(raw_pos_xls,1)-2)/(duration*60)
time_lapse=1000/FPS

num_of_tracks = (size(raw_pos_xls,2)-1)/3


%% saves the x-y coordinates in matrix form
raw_pos(1,1)=0; %sets 1st time point to 0
raw_pos(2:size(raw_pos_xls,1)-2,1)=(time_lapse)*(1:size(raw_pos_xls,1)-3); % fills in time

% for count_track=1:num_of_tracks
%     for count_frame=1:size(raw_pos_xls,1)-2
%         raw_pos(count_frame,count_track*2)=str2num(raw_pos_xls{count_frame+2,(count_track*3)-1}); %fills the x coordinates 
%         raw_pos(count_frame,(count_track*2)+1)=str2num(raw_pos_xls{count_frame+2,count_track*3}); %fills the y coordinates 
%     end
%     indi_raw_pos{1,count_track}=raw_pos(1:end,count_track*2:(count_track*2)+1);  %stores the x-y coordinates of each track separately
% end


for count_track=1:num_of_tracks
    for count_frame=1:size(raw_pos_xls,1)-2
        raw_pos(count_frame,count_track*2)=(raw_pos_xls{count_frame+2,(count_track*3)-1}); %fills the x coordinates 
        raw_pos(count_frame,(count_track*2)+1)=(raw_pos_xls{count_frame+2,count_track*3}); %fills the y coordinates 
    end
    indi_raw_pos{1,count_track}=raw_pos(1:end,count_track*2:(count_track*2)+1);  %stores the x-y coordinates of each track separately
end

    
corr_vs_win_size=[];
result=[];

for win_size = 2:size(raw_pos,1)
    win_size
    corr_vs_win_size(size(corr_vs_win_size,1)+1,1)=win_size/FPS;
    
    corr_matrices={};
    avg_corr_matrix(1:num_of_tracks,1:num_of_tracks)=0;
    for win_choice=1:size(raw_pos,1)-(win_size-1)
        
        parfor count_track=1:num_of_tracks
            indi_pos_0_mean{1,count_track}=indi_raw_pos{1,count_track}(win_choice:win_choice+win_size-1,:);
            mean_pos = mean (indi_pos_0_mean{1,count_track});
            indi_pos_0_mean{1,count_track}(:,1)=indi_pos_0_mean{1,count_track}(:,1)-mean_pos(:,1); %sets the mean X to 0
            indi_pos_0_mean{1,count_track}(:,2)=indi_pos_0_mean{1,count_track}(:,2)-mean_pos(:,2); %sets the mean Y to 0
        end
        
        
        for counter1=1:num_of_tracks
            for counter2=1:num_of_tracks
                correlation=corrcoef(indi_pos_0_mean{1,counter1},indi_pos_0_mean{1,counter2});
                corr_matrices{1,win_choice}(counter1,counter2)=correlation(1,2);
            end
        end
        
        avg_corr_matrix=avg_corr_matrix+corr_matrices{1,win_choice};
        
    end
    
    avg_corr_matrix = avg_corr_matrix./win_choice;
    
    linear_result=[];
    parfor count_track=1:num_of_tracks
        linear_result=horzcat(linear_result,avg_corr_matrix(count_track,count_track+1:end));
    end
    result=vertcat(result,linear_result);
    
    
end

corr_vs_win_size=horzcat(corr_vs_win_size,result);





positive_corr=result>=0;
n_positive=sum(positive_corr,2);
positive_corr=positive_corr.*result;
avg_pos_corr = sum(positive_corr,2)./n_positive;

negative_corr=result<=0;
n_negative=sum(negative_corr,2);
negative_corr=negative_corr.*result;
avg_neg_corr = sum(negative_corr,2)./n_negative;

figure
plot(corr_vs_win_size(:,1),avg_pos_corr)
hold on
plot(corr_vs_win_size(:,1),avg_neg_corr)
hold off