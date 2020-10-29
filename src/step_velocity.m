%% Clean start
clc; clear all; close all;

%% Initialization
startup;

%% Get data from file - sensorLog_4downwait - sensorLog_4upwait - sensorLog_7sdownwait - sensorLog_7supwait-mixl - sensorLog_fewstepsfull
filename = '../collected-data/sensorLog_fewstepsfull.txt';
usestreamlog=0;
showfigure=1;
showTable=1;
debugging=0;    
SamplePeriod=1/100; %period of time 10 ms to s /1000

%% Read data stream
import('com.liu.sensordata.*');  % Used to receive data.
reader = FileSensorDataReader(filename);
reader.start();
i = 1;

if usestreamlog == 1
    try
        %% Create data link
        server = StreamSensorDataReader(3400);
        % Makes sure to resources are returned.
        sentinel = onCleanup(@() server.stop());

        server.start();  % Start data reception.
        catch e
        fprintf(['Unsuccessful connecting to client!\n' ...
          'Make sure to start streaming from the phone *after* '...
                 'running this function!']);
        return;
    end
end

while reader.status()
    data = reader.getNext(5);  
    if ~isnan(data(1)) && ~any(isnan(data(2:10)))
        if i == 1
            time(i) = data(1,1);
        end
        if i > 1 && data(1,1) > time(i-1)
            time(i) = data(1,1);
            acc(:,i) = data(1,2:4);     %m/s^2
            gyr(:,i) = data(1,5:7);     %rad/s
            mag(:,i) = data(1,8:10);    %µT
            quat(:,i) = data(1,18:21);
            i = i + 1;
        end
        if i == 1
            i = i + 1;
        end
    end
end
t0 = time(1);
time = time - t0;
time = time/1000; % Convert from ms to s

% Plot the starting and ending vertical line
steps_range;

%% Plot RAW data   
if showfigure == 1
    figure;
    subplot(411);
    plot(time,acc(1,:),'b',time,acc(2,:),'g',time,acc(3,:),'r');
    xline(r_start);  xline(r_end);  
    title(['Accelerometer']);
    xlabel('Time [s]');
    ylabel('Acceleration [m/s^2]');
    legend('x','y','z');

    subplot(412);
    plot(time,gyr(1,:),'b',time,gyr(2,:),'g',time,gyr(3,:),'r');
    xline(r_start);  xline(r_end);  
    title(['Gyroscope']);
    xlabel('Time [s]');
    ylabel('Angular velocity [rad/s]');
    legend('x','y','z');

    subplot(413);
    plot(time,mag(1,:),'b',time,mag(2,:),'g',time,mag(3,:),'r');
    xline(r_start); xline(r_end);
    title(['Magnetometer']);
    xlabel('Time [s]');
    ylabel('Magnetic field [\mu T]');
    legend('x','y','z');

    subplot(414);
    plot(time,quat(1,:),'b',time,quat(2,:),'g',time,quat(3,:),'r',time,quat(4,:),'m');
    xline(r_start);  xline(r_end);  
    title( ['Orientation quaternion']);
    xlabel('Time [s]');
    ylabel('Orientation quaternion');
    legend('q0','q1','q2','q3');
end

%% Plot Modified data

% Detection method (Second one)
rms_steps_d;

% Estimated Position
pose_estimation;

% Classify the user behavior
%classify_behavior;

end