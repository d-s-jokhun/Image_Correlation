clear all

[~,~,raw_pos_xls]=xlsread('drift corrected_cir_001.xls');
duration=2 %in min
win_size_frames = 80 %in frames


FPS=(size(raw_pos_xls,1)-2)/(duration*60)
time_lapse=1/FPS; %in sec
win_size_time=time_lapse*win_size_frames

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



corr_matrices={};
avg_corr_matrix(1:num_of_tracks,1:num_of_tracks)=0;
for win_choice=1:size(raw_pos,1)-(win_size_frames-1)
    
    indi_pos_0_mean={};
    
    for count_track=1:num_of_tracks
        ean{1,count_track}=indi_raw_pos{1,count_track}(win_choice:win_choice+win_size_frames-1,:);
        indi_pos_0_mean{1,count_track}=indi_raw_pos{1,count_track}(win_choice:win_choice+win_size_frames-1,:);
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

figure
imagesc(avg_corr_matrix)
caxis([-1 1])
colormap jet
colorbar


%%get rid of all the autocorrelation (coeff=1) and find pair with highest coeff
mat4max_corr = avg_corr_matrix <1;
mat4max_corr = mat4max_corr.*avg_corr_matrix;
[maxi1,Idx1] = max(mat4max_corr);
[maxi2,Idx2] = max(maxi1);
max_corr_pair=[Idx2,Idx1(Idx2)]
max_corr_coeff = avg_corr_matrix(Idx2,Idx1(Idx2))


figure
for count = 1:num_of_tracks
plot (raw_pos(1:80,count*2),raw_pos(1:80,(count*2)+1),'k')
hold on
end
for count = [max_corr_pair]
plot (raw_pos(1:80,count*2),raw_pos(1:80,(count*2)+1),'r')
hold on
end
hold off



