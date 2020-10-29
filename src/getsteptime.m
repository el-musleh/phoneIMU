% Get the matrix value at the peak time detected
% This file suppose to change the time from s to time step T
% input= is the time like 12.65 s  (peak value)
% output= is the time(value) to get the same time 452 e.g time(452)=12.65

function out = getsteptime(input1, input2)
    for i=1:length(input1)
        if (input2+0.05 > input1(i)) && (input2-0.05 < input1(i))
            out=i;
            break;
        end
    end
end
