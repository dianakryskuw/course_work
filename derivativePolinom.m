function [ dpoli ] = derivativePolinom( pow, x, index )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


if(pow==1)
        m=size(x);
        for i=1:m(1)
            dpoli(i,1)=0;
        end;
        for i=2:m(2)+1
            if(i-1==index)
                dpoli(:,i)=1;
            else
                dpoli(:,i)=0;
            end;
        end;
else
        dpoli=derivativePolinom(pow-1,x,index);
        mp=size(dpoli);
        npol=mp(2);
        power=ceil(pow/2);
        j=npol+1;
        m=size(x);
        for i=1:m(2)
        for l=(i+1):m(2)
            for k=1:power
                if(i==index)
                    dpoli(:,j)=(pow-k).*x(:,i).^(pow-k-1).*x(:,l).^k;
                elseif(l==index)
                    dpoli(:,j)=k.*x(:,i).^(pow-k).*x(:,l).^(k-1);
                else
                    dpoli(:,j)=0;
                end;
                    j=j+1;
            end;
        end;
        end;
        for i=1:m(2)
            if(i==index)
                dpoli(:,j)=pow.*x(:,i).^(pow-1);
            else
                dpoli(:,j)=0;
            end;
                j=j+1;
        end;

end

