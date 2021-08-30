%%load input signal y y y y in a matrix called input_data



duration=2 %in min

FPS=size(input_data,1)/(duration*60)
time_lapse=1/FPS;


for count=1:size(input_data,2)
    autocorr_output(:,count)=autocorr(input_data(:,count),size(input_data,1)-1);
end

for count_col=1:size(autocorr_output,2)
    for count_row=1:size(autocorr_output,1)
        if autocorr_output(count_row,count_col)<=0
            decorr_time(count_col,1)=(count_row-1)*time_lapse;
            break
        end
    end
end

