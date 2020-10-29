%% Get steps from RMS signal (norm of acc)

% Control
get_checking_range=0.8;    %0.8    0.95    1.0 % time after the peak. to avoid wrong signals
get_bwpeak_range=0.8; % effect the number of peaks detected under the threshold!
chk_4end_d=0.3;  %0.51(20)  0.3(35+)

% Assign account norm value to each time step, not only one value
acc_norm_v=0;
for i=1:numel(time)
    acc_norm_v(i)=norm(acc(:,i));
end
plotfromzero=mean(acc_norm_v);
for i=1:numel(time)
    acc_norm_v(i)=acc_norm_v(i)-plotfromzero;
end

% Get all peaks from norm 3-axis
get_old_value=0;    get_old_time=0;
get_all_peaks=0;    get_counter=1;
get_all_time=0;
for i=r_s_c:r_e_c
    if abs(acc_norm_v(i)) > abs(get_old_value)
        get_old_value= acc_norm_v(i);
        get_old_time= time(i);
    elseif abs(acc_norm_v(i)) < abs(get_old_value)
       if  get_old_time+0.5 < time(i)
           get_all_peaks(get_counter)= get_old_value;
           get_all_time(get_counter)= get_old_time;
%            stem(get_old_time,get_old_value);
           get_counter=get_counter+1;
           get_old_value=0;
           get_old_time=0;
       end
    end
end

% Get only peaks above the threshold
get_right_peaks=0;  get_right_time=0; 
get_threshold_value=mean(get_all_peaks)*get_bwpeak_range;
get_counter=1;
for i=1:numel(get_all_peaks)
    if (abs(get_threshold_value) < abs(get_all_peaks(i)))
       get_right_peaks(get_counter)= get_all_peaks(i);
       get_right_time(get_counter)= get_all_time(i);
       get_counter=get_counter+1;
    end
end

% Get peaks with width ~1sec only
step_value=0;  step_time=0; 
get_counter=1;  get_check_done=1;
manched=0;  
while (get_check_done <= numel(get_right_time)-1)
    for i=get_check_done:numel(get_right_time)-1
        if get_right_time(i)-chk_4end_d > time(r_s_c) & get_right_time(i)+chk_4end_d < time(r_e_c)
           %if next one far add it
           if  get_right_time(i)+get_checking_range < get_right_time(i+1)
               step_time(get_counter)= get_right_time(i);
               step_value(get_counter)= get_right_peaks(i);
               get_counter=get_counter+1;
           %else close to each other, take the higer peak value!
           elseif get_right_time(i)+get_checking_range >= get_right_time(i+1) 
               if get_right_peaks(i) > get_right_peaks(i+1)
                   step_time(get_counter)= get_right_time(i);
                   step_value(get_counter)= get_right_peaks(i);
                   get_counter=get_counter+1;
                   %skip nex i and continue
                   get_check_done = i+2;
                   break;
               end
           end
        end
    end
    %check the last peak since not reached!
    if i==numel(get_right_time)-1
        if get_right_time(i)+chk_4end_d < time(r_e_c)
            %close to each other, take the best match!
            if get_right_time(i+1)-get_checking_range < get_right_time(i)
                %check if the time already added in the new list!
                for ii=1:numel(step_time)
                    if  get_right_time(i) == step_time(ii)
                        manched =1;
                        if get_right_peaks(i+1) > step_value(ii)    %get_counter-1
                            step_time(ii)= get_right_time(i+1);
                            step_value(ii)= get_right_peaks(i+1);
                        end
                        break;
                    end
                end
                if ~manched
                    if get_right_peaks(i+1) > get_right_peaks(i)
                        step_time(get_counter)= get_right_time(i+1);
                        step_value(get_counter)= get_right_peaks(i+1);
                    elseif get_right_peaks(i+1) < get_right_peaks(i)
                        step_time(get_counter)= get_right_time(i);
                        step_value(get_counter)= get_right_peaks(i);
                    end
                end
            %it's far so add it
            elseif get_right_time(i+1)-get_checking_range > get_right_time(i) 
                step_time(get_counter)= get_right_time(i+1);
                step_value(get_counter)= get_right_peaks(i+1);
            end
            break;
        end
    end
end

% Start/End for each step
drawaline=-1; my_i=1;    startig_my_i=0; rang_my_i=0;
while my_i <= size(time,2)
    repeat=0;
    my_i_v=0;
    if my_i <= r_s_c || my_i >= r_e_c
        drawaline(my_i)=-1;
        my_i=my_i+1;
    else
        if acc_norm_v(my_i) > -0.5 && acc_norm_v(my_i) < +0.5
            drawaline(my_i)=-1;
            my_i=my_i+1;
        elseif acc_norm_v(my_i) < -0.5 || acc_norm_v(my_i) > +0.5
            startig_my_i=my_i;
            while ~repeat
                repeat=1;
                my_i_v=my_i;
                if time(my_i_v)+0.15 < time(r_e_c)
                    rang_my_i=getsteptime(time,(time(my_i_v)+0.15));
                    for my_i=my_i_v:rang_my_i
                        if acc_norm_v(my_i) < -0.5 || acc_norm_v(my_i) > +0.5
                            repeat=0;
                            my_i=my_i+1;
                            break;
                        end
                    end
                else
                    my_i=my_i+1;
                    break;
                end
            end
            if (time(my_i)-time(startig_my_i)) > 0.7
                for i=startig_my_i:my_i
                    drawaline(i)=-6;
                end
            else
                for i=startig_my_i:my_i
                    drawaline(i)=-1;
                end
            end
            startig_my_i=0;
        end
    end
end

% Find pre_valley & pos_valley for each step!
pre_valley_t=0; pre_valley_v=0; count_pre_vv=1;
pos_valley_t=0;     pos_valley_v=0; count_pos_vv=1;
for i=1:size(time,2)-1
    if drawaline(i)==-1 && drawaline(i+1)==-6
        pre_valley_v(count_pre_vv)= acc_norm_v(i+1);
        pre_valley_t(count_pre_vv)= time(i+1);
        count_pre_vv=count_pre_vv+1;
    end
    if drawaline(i)==-6 && drawaline(i+1)==-1
        pos_valley_v(count_pos_vv)= acc_norm_v(i+1);
        pos_valley_t(count_pos_vv)= time(i+1);
        count_pos_vv=count_pos_vv+1;
    end
end

% Re-sort step_time
count_my_step_p=1; crt_list=0; count_crt_list=1; newlist_t=0;   newlist_v=0; newlsi_step_value=0; onebyonecunt=1;
for i=1:length(pre_valley_v)
    for ii=1:length(step_time)
        if step_time(ii) > pre_valley_t(i) && step_time(ii) < pos_valley_t(i)
            crt_list(count_crt_list)=step_time(ii);
            newlsi_step_value(count_crt_list)=step_value(ii);
            count_crt_list=count_crt_list+1;
        end
    end
    for i4=1:length(crt_list)
        if max(newlsi_step_value(:))== newlsi_step_value(i4)
            newlist_t(onebyonecunt)= crt_list(i4);
            newlist_v(onebyonecunt)= newlsi_step_value(i4);
        end
    end
    onebyonecunt=onebyonecunt+1;
    count_crt_list=1;
    crt_list=0;
    newlsi_step_value=0;
end
step_time=newlist_t;
step_value=newlist_v;

% Plot full peaks as if QRS-complex plotting! in the same figure
if showfigure == 1
    figure;
    subplot(311);
    hold on;
    plot(time,acc_norm_v,'b');
    plot(pre_valley_t,pre_valley_v,'rs','MarkerFaceColor','g');
    plot(step_time,step_value,'rv','MarkerFaceColor','r');
    plot(pos_valley_t,pos_valley_v,'rs','MarkerFaceColor','black');
    plot(time,drawaline,'-k');
    xline(r_start);  xline(r_end); 
    title(['RMS of 3-Axes']);
    xlabel('Time [s]');
    ylabel('Magnitude of Acceleration');
    legend('norm Acceleration','pre valley','paek','pos valley');
    hold off;
end

% Position and Velocity from cumtrapz
vx_z1=0;    vx_z2=0;    vx_z3=0;
px_z=0;     px_x=0;     px_y=0;
counter=time(r_s_c);
catchit=1;
if showfigure == 1
    hold on;
    for i=1:length(pre_valley_t)
        L1=getsteptime(time, pre_valley_t(i));
        L2=getsteptime(time, pos_valley_t(i));
        vx_z1=[vx_z1, SamplePeriod*cumtrapz(acc(1,[L1:L2]))];
        vx_z2=[vx_z2, SamplePeriod*cumtrapz(acc(2,[L1:L2]))];
        vx_z3=[vx_z3, SamplePeriod*cumtrapz(acc(3,[L1:L2]))];
        counter=[counter, time(L1:L2)];
    end
    counter=[counter, time(r_e_c)];
    vx_z1=[vx_z1, 0];   vx_z2=[vx_z2, 0];   vx_z3=[vx_z3, 0];
    px_x= SamplePeriod*cumtrapz(vx_z1);
    px_y= SamplePeriod*cumtrapz(vx_z2);
    px_z= SamplePeriod*cumtrapz(vx_z3);
    subplot(312);
    plot(counter,vx_z3,'b');    %counter,vx_z1,'r',counter,vx_z2,'g',counter,vx_z3,'b'
    xline(r_start);  xline(r_end);
    title('Velocity');
    xlabel('Time [s]');
    ylabel('Estimated Speed [m/s]');
    legend('z');    %'x','y','z'
    subplot(313);
    plot(counter,px_z,'b'); %counter,px_x,'r',counter,px_y,'g',counter,px_z,'b'
    title('Position'); 
    ylabel('Position Coordination [m]');
    xlabel('Time [s]');
    xline(r_start);  xline(r_end); 
    legend('z');    %'x','y','z'
    hold off;
end

% Get the max speed between each step
cpeakv=0;   cpeakt=0;   oldcpeakv=0; ncpv=0;    
for i=1:length(pre_valley_t)
    ncpv=ncpv+1;
    oldcpeakv=0;
    for ii=getsteptime(counter,pre_valley_t(i)):getsteptime(counter,pos_valley_t(i))
        if abs(oldcpeakv) < abs(vx_z3(ii))
            oldcpeakv=abs(vx_z3(ii));
            cpeakv(ncpv)=abs(oldcpeakv);
            cpeakt(ncpv)=counter(ii);
        end
    end
end

% Get the direction
counteritup=0;  counteritdown=0;    invpos=0;
for i=1:length(px_z)
    if px_z(i) > 0
        counteritup=counteritup+1;
    elseif px_z(i) < 0
        counteritdown=counteritdown+1;
    end
end
if counteritdown < counteritup
    invpos=1;
end