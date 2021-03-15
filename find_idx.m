function [midpoint_time_idx, time_integral] = find_idx(time, current, time_treshold)

time_temp = zeros(size(current));
time_integral = zeros(size(current));

idx = -0.1 < current & current < 0.1;
time_temp(idx) = time(idx)*60; % time in minutes
time_temp = [time_temp;0];     % Force time to zero to complete integration by adding an extra 0
 
int_temp = 0;
jj = 0;
for ii= 2: length(time_temp)
    if time_temp(ii)~= 0
        jj = jj+1;
        int_temp = int_temp+(time_temp(ii)-time_temp(ii-1));
    else
        time_integral(ii-jj:ii-1) = int_temp-time_temp(ii-jj);
        jj = 0;
        int_temp = 0;
    end
end

seg = time_integral>=time_treshold;
seg = [seg;0]; % Force segment to zero to complete integration by adding an extra 0
jj = 0;tt = 0;
for ii = 1: length(seg)
    if seg(ii)~=0
        jj = jj+1;
        if jj == 1
            tt = tt+1;
        end
    else
        if jj>0
            midpoint_time_idx(tt) = ii - round(jj/2);   
        end      
        jj = 0;
    end
end
%time_mid_point=time(midpoint_time_idx);
end

