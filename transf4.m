function x = transf( xp1,xn1,xp2,xn2,xp3,xn3,xp4,xn4, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x=lhsdesign(4,n);
x=x';
a1=xp1;
b1=xn1-xp1;
a2=xp2;
b2=xn2-xp2;
a3=xp3;
b3=xn3-xp3;
a4=xp4;
b4=xn4-xp4;
for i=1:length(x)
    x(i,1)=a1+b1*x(i,1);
    x(i,2)=a2+b2*x(i,2);
    x(i,3)=a3+b3*x(i,3);
    x(i,4)=a4+b4*x(i,4);
end;

end


