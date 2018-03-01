function x = transf( xp,xn, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

count=length(xp);
x=lhsdesign(count,n);
x=x';
for i=1:count
a(i)=xp(i);
b(i)=xn(i)-xp(i);
end;
for i=1:length(x)
for j=1:count
    x(i,j)=a(j)+b(j)*x(i,j);
end;
end;
end
