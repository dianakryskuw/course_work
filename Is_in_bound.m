function [ res ] = Is_in_bound( lower_bound, upper_bound, point )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
res=1;
for i=1:length(point)
    if ((point(i)<lower_bound(i))||(point(i)>upper_bound(i)))
        res=0;
        break;
    end;
end;
end

