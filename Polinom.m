function [ poli ] = Polinom( pow, x )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if(pow==1)
        m=size(x);
        for i=1:m(1)
            poli(i,1)=1;
        end;
        for i=2:m(2)+1
            poli(:,i)=x(:,i-1);
        end;
else
        poli=Polinom(pow-1,x);
        mp=size(poli);
        npol=mp(2);
        power=ceil(pow/2);
        j=npol+1;
        m=size(x);
        for i=1:m(2)
        for l=(i+1):m(2)
            for k=1:power
                poli(:,j)=x(:,i).^(pow-k).*x(:,l).^k;
                j=j+1;
            end;
        end;
        end;
        for i=1:m(2)
            poli(:,j)=x(:,i).^pow;
            j=j+1;
        end;
end

