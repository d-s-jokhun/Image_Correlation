%%% Written by D.S.Jokhun on 05/01/2017
%%% Finds the 2D or 3D cross correlation between different trajectories
%%% Finds cross correlation as a function of observational window
%%% Do not forget to change lines 10 to 13


clc
clear all

smallest_win_size =60; % in sec (5)
win_size_incrememt =5; % in sec (5)
tot_duration=2 %in min
filenames = dir (['Corrected_rec_*.xls']);


pre_XCorr_vs_ObsWin_TeloVar=[];
pre_XCorr_vs_ObsWin_TimeVar=[];
pre_AvgXCorr_vs_ObsWin=[];


for f=1:size(filenames,1)
    filenames(f).name
    
    
    [~,~,telo_pos_xls]=xlsread(filenames(f).name);
    raw_telo_pos = cell2mat (telo_pos_xls);
    
    FPS=(size(raw_telo_pos,1)-2)/(tot_duration*60);
    time_lapse=1000/FPS;  % in ms
    
    smallest_win_size_frames = round(smallest_win_size * FPS);  % smallest window size in terms of number of frames
    win_size_incrememt_frames = round (win_size_incrememt * FPS); % increment in window size in terms of number of frames 
    
    num_of_Tracks = (size(raw_telo_pos,2)-1)/3
    num_of_TimePoints = size(raw_telo_pos,1)-2;
    
    WS=0;
    for win_size=smallest_win_size_frames:win_size_incrememt_frames:num_of_TimePoints   % setting the smallest and the largest time window as well as the jump
        win_size/FPS
        WS=WS+1;
        XCorr_SegmentPairs_cell_n=[];
        
        for win_choice=1:num_of_TimePoints-(win_size-1)
            
            
            segment_pos=[];
            segment_pos (1:win_size,1:num_of_Tracks*3) = raw_telo_pos(3+(win_choice-1):3+(win_choice-1)+(win_size-1),2:end);  % X,Y,Z positions of all tracks from a specific window choice with a specific window size
            segment_mean_pos = mean(segment_pos);
    
    
            segment_0_mean=[];
            for column_count = 1:num_of_Tracks*3
                segment_0_mean(:,column_count)= segment_pos(:,column_count)- segment_mean_pos(1,column_count);
            end
            
            
            pair_num=0;
            for count1_TrackNum = 1:num_of_Tracks
                for count2_TrackNum = count1_TrackNum+1:num_of_Tracks
                    pair_num=pair_num+1;
                    corr_coff=corrcoef(segment_0_mean(:,(count1_TrackNum*3)-2:count1_TrackNum*3),segment_0_mean(:,(count2_TrackNum*3)-2:count2_TrackNum*3));
                    XCorr_SegmentPairs_cell_n (pair_num,win_choice)= corr_coff(1,2);    % each column corresponds to correlations between all the possible pairs for each window choice of a particular window size
                end
            end
        end
        
        pre_XCorr_vs_ObsWin_TeloVar(WS,(f*2)-1)=mean(mean(XCorr_SegmentPairs_cell_n));
        pre_XCorr_vs_ObsWin_TeloVar(WS,(f*2))= std(mean(XCorr_SegmentPairs_cell_n,2));    % SD reflecting how the correlation varies between different telomere pairs within that cell
        pre_XCorr_vs_ObsWin_TimeVar(WS,(f*2)-1)=mean(mean(XCorr_SegmentPairs_cell_n));
        pre_XCorr_vs_ObsWin_TimeVar(WS,(f*2))= std(mean(XCorr_SegmentPairs_cell_n,1));    % SD reflecting the repeatability in mean correlation of the cell at that window size by choosing different observation windows
        
        pre_AvgXCorr_vs_ObsWin(WS,f)=mean(mean(XCorr_SegmentPairs_cell_n));
        
        
        
    end
    

end


time_win(1:size(pre_XCorr_vs_ObsWin_TeloVar,1),1)=smallest_win_size+(((1:size(pre_XCorr_vs_ObsWin_TeloVar,1))-1)*win_size_incrememt);
result_XCorr_vs_ObsWin_TeloVar=horzcat(time_win,pre_XCorr_vs_ObsWin_TeloVar);
result_XCorr_vs_ObsWin_TimeVar=horzcat(time_win,pre_XCorr_vs_ObsWin_TimeVar);

result_AvgXCorr_vs_ObsWin(:,1)=time_win;
result_AvgXCorr_vs_ObsWin(:,2)=mean(pre_AvgXCorr_vs_ObsWin,2);
result_AvgXCorr_vs_ObsWin(:,3)=std(pre_AvgXCorr_vs_ObsWin,0,2);  % variation between different cells

result_Avg_intracell_SD_TeloVar(:,1)=time_win;
result_Avg_intracell_SD_TeloVar(:,2)=mean(pre_XCorr_vs_ObsWin_TeloVar(:,2:2:end),2);   % mean of intracell variance with observation window
result_Avg_intracell_SD_TeloVar(:,3)=std(pre_XCorr_vs_ObsWin_TeloVar(:,2:2:end),0,2);  % SD of the mean for intracell variance


