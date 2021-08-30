%%% Written by D.S.Jokhun on 04/01/2017
%%% Finds the 2D or 3D cross correlation between different trajectories


% clc
clear all
% clearvars -except result_tether_pos;

filenames = dir (['Corrected_cyto*.xls']);
sig_prob=0.001 % p-value below which telomere pairs will be considered to be significantly correlated (or anti-correlated)
very_posi_corr=0.25;  % correlation value above which the telomere pair will be taken up for timescale analyses (significantly very positively correlated)
smallest_win_size =2; % in sec
win_size_incrememt =0.5; % in sec
duration=120; %in sec

result_avg_CorrCoeff=[];
result_CorrCoeff=[];
result_Corr_bet_CorrCoef_N_dist=[];
result_relationship_bet_CorrCoef_N_dist=[];
result_num_of_sig_posi_corr = [];
result_AvgCosTheeta_vs_timescale={};  %timescale|<CosTheeta of TeloPair1>|<CosTheeta of TeloPair1>|<CosTheeta of TeloPair1>
result_AvgCosTheeta_vs_timescale_combined=[];   %timescale|<CosTheeta of TeloPair1>|<CosTheeta of TeloPair1>|<CosTheeta of TeloPair1>
result_AvgCosTheeta_vs_timescale_combined(:,1)=smallest_win_size:win_size_incrememt:duration;

for f=1:size(filenames,1)
    filenames(f).name
    
    position=[];
    position_0_mean=[];
    corr_mat=[];
    corr_sig_mat=[];
    rslt_corr_vs_dist=[];
    
    [~,~,telo_pos_xls]=xlsread(filenames(f).name);
    raw_telo_pos = cell2mat (telo_pos_xls);
    
    num_of_Tracks = (size(raw_telo_pos,2)-1)/3
    num_of_TimePoints = size(raw_telo_pos,1)-2;

    position (1:num_of_TimePoints,1:num_of_Tracks*3) = raw_telo_pos(3:end,2:end);  % X,Y,Z positions of all tracks at all time
    

    
   
    mean_pos = mean(position);
    
    
    for column_count = 1:num_of_Tracks*3
        position_0_mean(:,column_count)= position(:,column_count)- mean_pos(1,column_count);
    end
    
    corr_sig_map=[];
    corr_sig_map(1:num_of_Tracks,1:num_of_Tracks)=2;
    telo_pair=0;
    sig_posi_telo_pair=0;
    sig_very_posi_telo_pair=0;
    AvgCosTheeta_vs_timescale_f=[];
    AvgCosTheeta_vs_timescale_f(:,1)=smallest_win_size:win_size_incrememt:duration;
    for count1_TrackNum = 1:num_of_Tracks
        for count2_TrackNum = count1_TrackNum+1:num_of_Tracks
            telo_pair=telo_pair+1;
            [corr_coff,corr_sig]=corrcoef(position_0_mean(:,(count1_TrackNum*3)-2:count1_TrackNum*3),position_0_mean(:,(count2_TrackNum*3)-2:count2_TrackNum*3));
            corr_mat(count1_TrackNum,count2_TrackNum)=corr_coff(1,2);
            corr_sig_mat(count1_TrackNum,count2_TrackNum)=corr_sig(1,2);
            if corr_sig(1,2)<sig_prob
                corr_sig_map(count1_TrackNum,count2_TrackNum)=0;
                if corr_coff(1,2)> 0
                    sig_posi_telo_pair=sig_posi_telo_pair+1;
                end
                if corr_coff(1,2)> very_posi_corr
                    sig_very_posi_telo_pair=sig_very_posi_telo_pair+1;
                    % investigating the correlation between significantly very positive telomere pairs as a function of timescale
                    track_for_timescale_studies=[];
                    track_for_timescale_studies=[position_0_mean(:,(count1_TrackNum*3)-2:count1_TrackNum*3),position_0_mean(:,(count2_TrackNum*3)-2:count2_TrackNum*3)];  % XYZ coordinates of telomere 1(column 1-3), XYZ coordinates of telomere 2(column 4-6) with time (rows)
                    
                    time_lapse=duration/(size(track_for_timescale_studies,1)-1);
                    WS=0;
                    for WinSize=smallest_win_size:win_size_incrememt:duration
                        WS=WS+1;
                        possible_displacement_vectors=track_for_timescale_studies(round(WinSize/time_lapse)+1:end,:)-track_for_timescale_studies(1:end-round(WinSize/time_lapse),:);
                        % cos(theeta)= (A.B) / |A|*|B|
                        numerators=dot(possible_displacement_vectors(:,1:3),possible_displacement_vectors(:,4:6),2);
                        denominators=sqrt((sum((possible_displacement_vectors(:,1:3).^2),2)).*(sum((possible_displacement_vectors(:,4:6).^2),2)));
                        CosTheeta=numerators./denominators;
                        AvgCosTheeta_vs_timescale_f(WS,sig_very_posi_telo_pair+1)=mean(CosTheeta);

                    end
                end
            else
                corr_sig_map(count1_TrackNum,count2_TrackNum)=1;
            end
            
            result_CorrCoeff(end+1,1)=corr_coff(1,2);   %gives the correlation coefficients between all the telomere pairs studied
            result_CorrCoeff(end,2)=corr_sig(1,2);   % gives the corresponding p-value to check for significance (p<0.05) mean that the correlaiton is significant
            rslt_corr_vs_dist(size(rslt_corr_vs_dist,1)+1,1)= sqrt(sum((mean_pos((count2_TrackNum*3)-2:count2_TrackNum*3)-mean_pos((count1_TrackNum*3)-2:count1_TrackNum*3)).^2));
            rslt_corr_vs_dist(size(rslt_corr_vs_dist,1),2)=corr_coff(1,2);
        end
    end
    
    result_AvgCosTheeta_vs_timescale{f}=AvgCosTheeta_vs_timescale_f;
    result_AvgCosTheeta_vs_timescale_combined(1:size(AvgCosTheeta_vs_timescale_f,1),end+1:end+((size(AvgCosTheeta_vs_timescale_f,2))-1))=AvgCosTheeta_vs_timescale_f(:,2:end);
    
    
    
    result_num_of_sig_posi_corr(end+1,1:2)=[sig_posi_telo_pair,telo_pair];  % 1st column is the number of sig posi telo pairs and the 2nd column is the total num of telo pairs in that cell.
    
    for count_fill = 1:num_of_Tracks
        corr_mat(count_fill,count_fill)=1;
        corr_sig_mat(count_fill,count_fill)=1;
    end
    
    
    
    fit_corr_vs_dist = fit(rslt_corr_vs_dist(:,1),rslt_corr_vs_dist(:,2),'poly1','Robust','Bisquare'); 
    coeffvals_corr_vs_dist = coeffvalues(fit_corr_vs_dist); % gives the coefficients of the line y=Mx+C  (M and C)
    

%     figure
%     imagesc(corr_mat)
%     caxis([-1 1])
%     colormap jet
%     colorbar
%     ax=gca;
%     ax.YDir = 'normal';
%     
%     c=[0 1 0;0 0 0;1 1 1];
%     figure
%     imagesc(corr_sig_map)
%     caxis([0 2])
%     colormap (c)
%     colorbar
%     ax=gca;
%     ax.YDir = 'normal';
%     
%     
%     
%     figure
%     scatter(rslt_corr_vs_dist(:,1),rslt_corr_vs_dist(:,2))
%     hold on
%     plot(fit_corr_vs_dist)
%     hold off
%     
    
    
    
    result_avg_CorrCoeff(size(result_avg_CorrCoeff,1)+1,1)=mean(rslt_corr_vs_dist(:,2));
    pre_corr_coeff_corrNdist=corrcoef(rslt_corr_vs_dist(:,1),rslt_corr_vs_dist(:,2));
    result_Corr_bet_CorrCoef_N_dist(size(result_Corr_bet_CorrCoef_N_dist,1)+1,1) = pre_corr_coeff_corrNdist(1,2);
    result_relationship_bet_CorrCoef_N_dist(size(result_relationship_bet_CorrCoef_N_dist,1)+1,1) = coeffvals_corr_vs_dist(2);   % Y-intercept
    result_relationship_bet_CorrCoef_N_dist(size(result_relationship_bet_CorrCoef_N_dist,1),2) = coeffvals_corr_vs_dist(1);   % gradient
    result_relationship_bet_CorrCoef_N_dist(size(result_relationship_bet_CorrCoef_N_dist,1),3) = -coeffvals_corr_vs_dist(2)/coeffvals_corr_vs_dist(1);   % X-intercept

    
    
end



result_sig_posi_corr_val=[];
posi_corr=result_CorrCoeff>0;
sig_corr=result_CorrCoeff<sig_prob;
for count=1:size(result_CorrCoeff,1)
    if posi_corr(count,1)==1
        if sig_corr(count,2)==1
            sig_posi_corr(count,1:2)=1;
            result_sig_posi_corr_val(end+1,1:2)=result_CorrCoeff(count,1:2);
        end
    end
end
if size(sig_posi_corr,1)<size(result_CorrCoeff,1)
    sig_posi_corr(size(result_CorrCoeff,1),1:2)=0;
end


Overall_frac_of_sig_posi_corr = (sum(sig_posi_corr(:,1)))/size(result_CorrCoeff,1)

result_Mean_of_AvgCosTheeta_vs_timescale=[];
result_Mean_of_AvgCosTheeta_vs_timescale(:,1)=result_AvgCosTheeta_vs_timescale_combined(:,1);
result_Mean_of_AvgCosTheeta_vs_timescale(:,2)=mean(result_AvgCosTheeta_vs_timescale_combined(:,2:end),2);
result_Mean_of_AvgCosTheeta_vs_timescale(:,3)=std(result_AvgCosTheeta_vs_timescale_combined(:,2:end),0,2);
result_Mean_of_AvgCosTheeta_vs_timescale(:,4)=median(result_AvgCosTheeta_vs_timescale_combined(:,2:end),2);



%%% Timescale dependence of CosTheeta (Gradients between 15 and 45sec)
result_TimescaleDependence=[];
%% gradient of mean
x=[];
y=[];
x=result_Mean_of_AvgCosTheeta_vs_timescale(27:87,1);
y=result_Mean_of_AvgCosTheeta_vs_timescale(27:87,2);
% PCA
[coeff,~,~,~,~] = pca([x,y]);
result_TimescaleDependence(1,1)=(coeff(2,1)/coeff(1,1));
% linear fit
[fit_AvgCosTheeta_vs_timescale] = fit(x,y,'poly1','Robust','Bisquare'); 
coeffvals_AvgCosTheeta_vs_timescale = coeffvalues(fit_AvgCosTheeta_vs_timescale); % gives the coefficients of the line y=Mx+C  (M and C)
result_TimescaleDependence(2,1)=coeffvals_AvgCosTheeta_vs_timescale(1);
%% gradient of fit of combined
x=[];
y=[];
for count=1:size(result_AvgCosTheeta_vs_timescale_combined,2)-1;
    x=vertcat(x,result_AvgCosTheeta_vs_timescale_combined(27:87,1));
    y=vertcat(y,result_AvgCosTheeta_vs_timescale_combined(27:87,count+1));
end
% PCA
[coeff,~,~,~,~] = pca([x,y]);
result_TimescaleDependence(3,1)=(coeff(2,1)/coeff(1,1));
% linear fit
[fit_AvgCosTheeta_vs_timescale] = fit(x,y,'poly1','Robust','Bisquare'); 
coeffvals_AvgCosTheeta_vs_timescale = coeffvalues(fit_AvgCosTheeta_vs_timescale); % gives the coefficients of the line y=Mx+C  (M and C)
result_TimescaleDependence(4,1)=coeffvals_AvgCosTheeta_vs_timescale(1);



















