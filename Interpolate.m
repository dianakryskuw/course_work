function [ xmin, fmin, count ] = Interpolate( f, str ,np,startPoint, lb, ub)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

narg=length(lb);
for i=1:narg
xl(i)=lb(i);
xu(i)=ub(i);
end

% % if(narg==2)
% %     Z1=xl(1):0.1:xu(1);
% %     Z2=xl(2):0.1:xu(2);
% %     [z1,z2]=meshgrid(Z1,Z2);
% %     mesh(z1,z2,f(z1,z2));
% % end;
% f=str2func(fun);


k=0.5;
    r=startPoint;
eps=0.001;
hold on
x1=[];
addx=[];
addy=[];
y=[];
n=np;
l=0;

for i=1:narg
coeff(i)=0.2*(xu(i)-xl(i));
end;
fInvokeCount=0;


while ((l==0)||(norm(rn-r)>1)||(abs(phio(r)-f(r))>eps))
    if(l>0)
        for i=1:narg
            coeff(i)=abs(rn(i)-r(i))*1.5;
        end;
    end;
      rn=r;

      
    for i=1:length(r)
lower_bound(i)=max(r(i)-coeff(i),xl(i));
upper_bound(i)=min(r(i)+coeff(i),xu(i));
    end;
      
%plot3(r(1),r(2),f(r(1),r(2)),'*r');
xx=transf(lower_bound,upper_bound,n);
% 
% x1d=x(:,1);
% x2d=x(:,2);
% x3d=x(:,3);

m=size(xx);

for i=1:m(1)
        myx=[xx(i,:)];
    yy(i)=f(myx);
    fInvokeCount=fInvokeCount+1;
end;

if(l>0)
j=1;
addx=[];
addy=[];
for i=1:m(1)
    if(Is_in_bound(lower_bound, upper_bound, x(i,:)))
        for k=1:narg
            addx(j,k)=x(i,k);
        end;
            addy(j)=y(i);
            j=j+1;
    end;
end;
end;

if(~isempty(addx))
x=[xx;addx];
else
    x=xx;
end;
if(~isempty(addy))
y=[yy';addy'];
else
    y=yy';
end;


if (strcmp(str,'lin')==1)
    X=Polinom(1,x);
   % X=[ x(:,1).^0, x(:,1), x(:,2)];
elseif (strcmp(str,'quad')==1)
    X=Polinom(2,x)
   % XX=[ x(:,1).^0,  x(:,1), x(:,2), x(:,1).*x(:,2), x(:,1).^2, x(:,2).^2];
elseif (strcmp(str,'cub')==1)
    X=Polinom(3,x);
    %X=[ x1.^0, x1, x2, x1.*x2, x1.^2, x2.^2, x2.*x1.^2, x1.*x2.^2, x1.^3, x2.^3];
elseif (strcmp(str,'quar')==1)
    X=Polinom(4,x);
   % X=[ x1.^0, x1, x2, x1.*x2, x1.^2, x2.^2, x1.*x2.^2, x1.*x2.^2, x1.^3, x2.^3, x1.*x2.^3, x2.*x1.^3, x1.^2.*x2.^2, x1.^4, x2.^4];
end;


Xt=X';
A=Xt*X;
b=Xt*y;
a=A\b;


if (strcmp(str,'lin')==1)
    phio=@(x)a.*Polinom(1,x);
    phio1=@(x)sum(Meta_model(x,a,1));
    phi=@(x1,x2)a(1)+a(2).*x1+a(3).*x2;
    phio2=@(x)a(1)+a(2).*x(1)+a(3).*x(2);
elseif (strcmp(str,'quad')==1)
    phio=@(x)sum(a'.*Polinom(2,x));
    phio1=@(x)sum(Meta_model(x,a,2));
    phi=@(x1,x2)a(1)+a(2).*x1+a(3).*x2+a(4).*x1.*x2+a(5).*x1.^2+a(6).*x2.^2;
    phio2=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(1).*x(2)+a(5).*x(1).^2+a(6).*x(2).^2;
elseif (strcmp(str,'cub')==1)
    phio=@(x)sum(a'.*Polinom(3,x));
    phio1=@(x)sum(Meta_model(x,a,3));
    phi=@(x1,x2)a(1)+a(2).*x1+a(3).*x2+a(4).*x1.*x2+a(5).*x1.^2+a(6).*x2.^2+a(7).*x2.*x1.^2+a(8).*x1.*x2.^2+a(9).*x1.^3+a(10).*x2.^3;
    phio2=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(1).*x(2)+a(5).*x(1).^2+a(6).*x(2).^2+a(7).*x(2).*x(1).^2+a(8).*x(1).*x(2).^2+a(9).*x(1).^3+a(10).*x(2).^3;
elseif (strcmp(str,'quar')==1)
    phio=@(x)Polinom(4,x).*a;
    phio1=@(x)sum(Meta_model(x,a,4));
    phi=@(x1,x2)a(1)+a(2).*x1+a(3).*x2+a(4).*x1.*x2+a(5).*x1.^2+a(6).*x2.^2+a(7).*x2.*x1.^2+a(8).*x1.*x2.^2+a(9).*x1.^3+a(10).*x2.^3+a(11).*x1.*x2.^3+a(12).*x2.*x1.^3+a(13).*x1.^2.*x2.^2+a(14).*x1.^4+a(15).*x2.^4;
    phio2=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(1).*x(2)+a(5).*x(1).^2+a(6).*x(2).^2+a(7).*x(2).*x(1).^2+a(8).*x(1).*x(2).^2+a(9).*x(1).^3+a(10).*x(2).^3+a(11).*x(1).*x(2).^3+a(12).*x(2).*x(1).^3+a(13).*x(1).^2.*x(2).^2+a(14).*x(1).^4+a(15).*x(2).^4;
end;
yy=[1,1];
phio(yy)
phio2(yy)
if(abs(upper_bound(1)-lower_bound(1))<eps)||(abs(coeff(1))<eps)
    break;
end;
% % % if(a1>a2)
% % % ta=a1;
% % % a1=a2;
% % % a2=ta;
% % % k=k+0.02;
% % % end;
% % % 
% % % if(b1>b2)
% % % tb=b1;
% % % b1=b2;
% % % b2=tb;
% % % k=k+0.02;
% % % end;
[r,out,opt]=fmincon(phio,rn,[],[],[],[],lower_bound,upper_bound)



% % % if(abs(phio(r)-My(r))<1)
% % % k=k-0.02;
% % % eps=0.00001;
% % % end;
% % % % 
% % % % if(r(1)<x1l)
% % % %     r(1)=x1l;
% % % % end;
% % % % if(r(1)>x1u)
% % % %     r(1)=x1u;
% % % % end;
% % % % if(r(2)<x2l)
% % % %     r(2)=x2l;
% % % % end;
% % % % if(r(2)>x2u)
% % % %     r(2)=x2u;
% % % % end;
l=l+1;
end;
disp('');
disp('');
disp('------------------------------RESULTS---------------------------------------------------------------');
disp('');
disp('');
xmin=r
fmin=f(xmin)
count=l
%plot3(r(1),r(2),f(r(1),r(2)),'ok');

%mesh(z1,z2,phi(z1,z2));
%hold off

end

