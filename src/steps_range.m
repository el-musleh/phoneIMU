%% Get steps detection region

% To ignore the values at the start and end, it's noise and to avoid count 
% it as step, this is only working for the log file.
r_s_c=0;            r_e_c=0;
r_start=0;          r_end=0;
counter_abc=0;      exit_2nd_loop=0;
notdirectly=0;      recent_mag_v=0;
pchange_end=0;

% Change the detection range according to the time lenght
TotalSampleTime= time(numel(time));
pre_eval = round(-0.5*TotalSampleTime+30);
pre_range_end= round(+3.5*TotalSampleTime+5);
pos_range_end= round(2.75*TotalSampleTime+32);
pos_eval = round(1.25*TotalSampleTime+46);

e_st_t= round((numel(time)*pre_eval)/100);
ve_ed_t = round((numel(time)*pos_eval)/100);
rangelong1=12;
rangelong2 =4;

% For start range vertical line
% from norm acc get the start peak!
oldpeak_time=0; oldpeak_value=0;
r_s_c_f=0;  r_start_f=0;
for i=1:e_st_t
    if norm(acc(:,i)) >= oldpeak_value
        oldpeak_value = norm(acc(:,i));
        oldpeak_time=i;
    end
end
r_s_c_f = oldpeak_time+pre_range_end;   
r_start_f = time(r_s_c_f);
% from mag value
found=0;
for i=e_st_t:-1:1
    if mag(2,i) > 0 & r_start == 0
        pchange_end = mag(2,i);
        recent_mag_v=mag(2,i);
        for ii=i:-1:1
            if pchange_end > mag(2,ii) & recent_mag_v >= mag(2,ii) & ~found & ~exit_2nd_loop
                counter_abc = counter_abc +1;
                if(mag(2,ii) < 0) & counter_abc > rangelong1 & ~found & ~exit_2nd_loop
                    recent_mag_v=mag(2,ii);
                    for iii=ii:-1:1
                        if recent_mag_v >= mag(2,iii)
                            notdirectly = notdirectly+1;
                            if notdirectly > 4
                                r_s_c = i+pre_range_end;
%                                 sprintf('y')
                                r_start = time(r_s_c);  
                                found = 1;
                            end
                        elseif(recent_mag_v < mag(2,iii))
                            notdirectly=0;
                            exit_2nd_loop = 1;
                            break;
                        end
                        recent_mag_v=mag(2,iii);
                    end
                end
            elseif recent_mag_v < mag(2,ii) | found | exit_2nd_loop
%                 sprintf(num2str(counter_abc))
                counter_abc=0;
                break;
            end
            recent_mag_v=mag(2,ii);
        end
        exit_2nd_loop=0;
    end
    if mag(1,i) > 0 & r_start == 0
        pchange_end = mag(1,i);
        recent_mag_v=mag(1,i);
        for ii=i:-1:1
            if pchange_end > mag(1,ii) & recent_mag_v >= mag(1,ii) & ~found & ~exit_2nd_loop
                counter_abc = counter_abc +1;
                if(mag(1,ii) < 0) & counter_abc > rangelong1 & ~found & ~exit_2nd_loop
                    recent_mag_v=mag(1,ii);
                    for iii=ii:-1:1
                        if recent_mag_v >= mag(1,iii)
                            notdirectly = notdirectly+1;
                            if notdirectly > 4
                                r_s_c = i+pre_range_end;
%                                 sprintf('x')
                                r_start = time(r_s_c);
                                found = 1;
                            end
                        elseif(recent_mag_v < mag(1,iii))
                            notdirectly=0;
                            exit_2nd_loop = 1;
                            break;
                        end
                        recent_mag_v=mag(1,iii);
                    end
                end
            elseif recent_mag_v < mag(1,ii) | found | exit_2nd_loop
%                 sprintf(num2str(counter_abc))
                counter_abc=0;
                break;
            end
            recent_mag_v=mag(1,ii);
        end
        exit_2nd_loop=0;
    end
    if mag(3,i) > 0 & r_start == 0
        pchange_end = mag(3,i);
        recent_mag_v=mag(3,i);
        for ii=i:-1:1
            if pchange_end > mag(3,ii) & recent_mag_v >= mag(3,ii) & ~found & ~exit_2nd_loop
                counter_abc = counter_abc +1;
                if(mag(3,ii) < 0) & counter_abc > rangelong1 & ~found & ~exit_2nd_loop
                    recent_mag_v=mag(3,ii);
                    for iii=ii:-1:1
                        if recent_mag_v >= mag(3,iii)
                            notdirectly = notdirectly+1;
                            if notdirectly > 4
                                r_s_c = i+pre_range_end;
%                                 sprintf('z')
                                r_start = time(r_s_c);
                                found = 1;
                            end
                        elseif(recent_mag_v < mag(3,iii))
                            notdirectly=0;
                            exit_2nd_loop = 1;
                            break;
                        end
                        recent_mag_v=mag(3,iii);
                    end
                end
            elseif recent_mag_v < mag(3,ii) | found | exit_2nd_loop
%                 sprintf(num2str(counter_abc))
                counter_abc=0;                
                break;
            end
            recent_mag_v=mag(3,ii);
        end
        exit_2nd_loop=0;
    end
end
% check the difference if big take the smaller  - check if anyone is zero
v_difference=abs(r_s_c_f-r_s_c);
% sprintf('r_s_c=%d  r_s_c_f=%d \n r_start=%d  r_start_f=%d\n %d',r_s_c, r_s_c_f, r_start, r_start_f, v_difference)
if v_difference > 50 & v_difference <= 100
    if r_s_c > r_s_c_f
        r_s_c=r_s_c_f+round(v_difference/2)-30;
        r_start= time(r_s_c);
    else
        r_s_c=r_s_c+round(v_difference/2)-30;
        r_start= time(r_s_c);
    end
elseif v_difference > 100
    if r_s_c > r_s_c_f
        r_s_c=r_s_c_f+round(v_difference/2)-80;
        r_start= time(r_s_c);
    else
        r_s_c=r_s_c+round(v_difference/2)-80;
        r_start= time(r_s_c);
    end
end
% sprintf('r_s_c=%d  r_start=%d ',r_s_c, r_start)
% else from % time
if r_start == 0
    for i=1:e_st_t
        if time(i) >= (2.1) & time(i) <= (2.4)
            r_s_c=i;
            r_start=time(i);
        end
    end
end

% For last range vertical line
found = 0;
while ~found & pos_eval > 49
    for i=numel(time):-1:ve_ed_t
        ve_ed_t = round((numel(time)*pos_eval)/100);
        if mag(2,i) < 0 & r_end == 0
            pchange_end = mag(2,i);
            recent_mag_v=mag(2,i);
            for ii=i:-1:ve_ed_t
                if pchange_end < mag(2,ii) & recent_mag_v <= mag(2,ii) & ~found & ~exit_2nd_loop
                    counter_abc = counter_abc +1;
                    if(mag(2,ii) > 0) & counter_abc > rangelong2 & ~found & ~exit_2nd_loop
                        recent_mag_v=mag(2,ii);
                        for iii=ii:-1:ve_ed_t
                            if recent_mag_v <= mag(2,iii)
                                notdirectly = notdirectly+1;
                                if notdirectly > 4
                                    r_e_c = iii-pos_range_end;
                                    r_end = time(r_e_c);
                                    found = 1;
                                end
                            elseif(recent_mag_v > mag(2,iii))
                                notdirectly=0;
                                exit_2nd_loop = 1;
                                break;
                            end
                            recent_mag_v=mag(2,iii);
                        end
                    end
                elseif recent_mag_v > mag(2,ii) | found | exit_2nd_loop
    %                 sprintf(num2str(counter_abc))
                    counter_abc=0;
                    break;
                end
                recent_mag_v=mag(2,ii);
            end
            exit_2nd_loop=0;
        end
        if mag(1,i) < 0 & r_end == 0
            pchange_end = mag(1,i);
            recent_mag_v=mag(1,i);
            for ii=i:-1:ve_ed_t
                if pchange_end < mag(1,ii) & recent_mag_v <= mag(1,ii) & ~found & ~exit_2nd_loop
                    counter_abc = counter_abc +1;
                    if(mag(1,ii) > 0) & counter_abc > rangelong2 & ~found & ~exit_2nd_loop
                        recent_mag_v=mag(1,ii);
                        for iii=ii:-1:ve_ed_t
                            if recent_mag_v <= mag(1,iii)
                                notdirectly = notdirectly+1;
                                if notdirectly > 4
                                    r_e_c = iii-pos_range_end; 
                                    r_end = time(r_e_c); 
                                    found = 1;
                                end
                            elseif(recent_mag_v > mag(1,iii))
                                notdirectly=0;
                                exit_2nd_loop = 1;
                                break;
                            end
                            recent_mag_v=mag(1,iii);
                        end
                    end
                elseif recent_mag_v > mag(1,ii) | found | exit_2nd_loop
    %                 sprintf(num2str(counter_abc))
                    counter_abc=0;
                    break;
                end
                recent_mag_v=mag(1,ii);
            end
            exit_2nd_loop=0;
        end
        if mag(3,i) < 0 & r_end == 0
            pchange_end = mag(3,i);
            recent_mag_v=mag(3,i);
            for ii=i:-1:ve_ed_t
                if pchange_end < mag(3,ii) & recent_mag_v <= mag(3,ii) & ~found & ~exit_2nd_loop
                    counter_abc = counter_abc +1;
                    if(mag(3,ii) > 0) & counter_abc > rangelong2 & ~found & ~exit_2nd_loop
                        recent_mag_v=mag(3,ii);
                        for iii=ii:-1:ve_ed_t
                            if recent_mag_v <= mag(3,iii)
                                notdirectly = notdirectly+1;
                                if notdirectly > 4
                                    r_e_c = iii-pos_range_end; 
                                    r_end = time(r_e_c);
                                    found = 1;
                                end
                            elseif(recent_mag_v > mag(3,iii))
                                notdirectly=0;
                                exit_2nd_loop = 1;
                                break;
                            end
                            recent_mag_v=mag(3,iii);
                        end
                    end
                elseif recent_mag_v > mag(3,ii) | found | exit_2nd_loop
    %                 sprintf(num2str(counter_abc))
                    counter_abc=0;
                    break;
                end
                recent_mag_v=mag(3,ii);
            end
            exit_2nd_loop=0;
        end
    end
    pos_eval=pos_eval-10;
end
if r_end == 0
    r_e_c = ve_ed_t;
    r_end=time(r_e_c);
end
