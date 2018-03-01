function [ arr ] = Meta_model( y,a,pow )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(pow==1)
        m=length(y)
            arr(1)=a(1);
        for i=2:m+1
            arr(i)=y(i-1).*a(i);
        end;
else
        arr=Meta_model(y,a,pow-1);
        npol=length(arr);
        power=ceil(pow/2);
        j=npol+1;
        m=length(y);
        for i=1:m
        for l=(i+1):m
            for k=1:power
                arr(j)=a(j).*y(i).^(pow-k).*y(l).^k;
                j=j+1;
            end;
        end;
        end;
        for i=1:m
            arr(j)=a(j).*y(i).^pow;
            j=j+1;
        end;
end

