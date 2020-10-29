%% Get Estimated Position

% Get all the steps time for plotting vertical line
all_vally_tt=0;
for i=1:length(pre_valley_t)
    if i>1
        all_vally_tt= [all_vally_tt, pre_valley_t(i), pos_valley_t(i)];
    else
        all_vally_tt= [pre_valley_t(i), pos_valley_t(i)];
    end
end

% Plot 2D figure show steps elevation
g=[0;0;9.81]; %Gravity of Earth m/s^2
a_est = zeros(3,1);
v_est = zeros(3,1);
p_est = zeros(3,1);
stime=0;
t=1;
for i=1:numel(pre_valley_t)
    for ii=(getsteptime(time,pre_valley_t(i))):(getsteptime(time,pos_valley_t(i)))
        q = [quat(1,ii);quat(2,ii);quat(3,ii);quat(4,ii)]';
        R = [2*(q(1)^2+q(2)^2)-1, 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
             2*(q(2)*q(3)+q(1)*q(4)), 2*(q(1)^2 +q(3)^2)-1, 2*(q(3)*q(4)-q(1)*q(2));
             2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 2*(q(1)^2+q(4)^2)-1];
        a_est(:,t) = (R)*acc(:,ii)-g;
        if isnan(a_est(:,t))        %if value NaN, use recent value
            a_est(:,t) = a_est(:,t-1);
        end
        v_est(:,t+1)=SamplePeriod*a_est(:,t)+v_est(:,t);
        p_est(:,t+1)=p_est(:,t)+SamplePeriod*v_est(:,t)+((SamplePeriod^2)/2)*a_est(:,t); 
        stime(t)=time(ii);
        t=t+1;
    end
end
stime(t)=time(r_e_c);

figure;
hold on;
subplot(311);
a_est(:,t) = (R.')*acc(:,ii+1)-g;
plot(stime,a_est(3,:),'b');   %stime,a_est(1,:),'r',stime,a_est(2,:),'g',stime,a_est(3,:),'b'
xline(r_start); xline(r_end); 
for i=1:length(all_vally_tt)
    xline(all_vally_tt(i));
end
title('Acceleration');    
xlabel('time [s]');
ylabel('Accelerometer [m/s^2]');
legend('z'); %'x','y','z'

subplot(312);
plot(stime,v_est(3,:),'b');     %stime,v_est(1,:),'r',stime,v_est(2,:),'r',stime,v_est(3,:),'b'
xline(r_start); xline(r_end); 
for i=1:length(all_vally_tt)
    xline(all_vally_tt(i));
end
title('Velocity');
xlabel('Time [s]');
ylabel('Estimated Speed [m/s]');
legend('z');    %'x','y','z'

subplot(313);
newposest=0;
if invpos==1
    for i=1:size(p_est,2)
        newposest(1, i)=abs(p_est(1,i));
        newposest(2, i)=abs(p_est(2,i));
        newposest(3, i)=abs(p_est(3,i));
    end
    plot(stime, newposest(3,:), 'b');    %stime, newposest(1,:),'r', stime,newposest(2,:), 'g',stime, newposest(3,:), 'b'
    xline(r_start); xline(r_end); 
    for i=1:length(all_vally_tt)
        xline(all_vally_tt(i));
    end
    title('Position');
    xlabel('Time [s]');
    ylabel('Position Coordination [m]');
    legend('z');    %'x','y','z'
    hold off;
else
    plot(stime, p_est(3,:), 'b');    %stime, p_est(1,:),'r', stime,p_est(2,:), 'g',stime, p_est(3,:), 'b'
    xline(r_start); xline(r_end); 
    for i=1:length(all_vally_tt)
        xline(all_vally_tt(i));
    end
    title('Position');
    xlabel('Time [s]');
    ylabel('Position Coordination [m]');
    legend('z');    %'x','y','z'
    hold off;
end
