function [ xmin, fmin, count ] = Interpolate( f, str ,np,startPoint, lb, ub,df1,df2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

narg=length(lb);
for i=1:narg
xl(i)=lb(i);
xu(i)=ub(i);
end;

k=0.5;
r=startPoint;
eps=0.00001;
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
      
    
xx=transf(lower_bound,upper_bound,n);

m=size(xx);

for i=1:m(1)
    myx=[xx(i,:)];
    yy(i)=f(myx);
    dfy1(i)=df1(myx);
    dfy2(i)=df2(myx);
    fInvokeCount=fInvokeCount+1;
end;

%yy=[yy,dfy1,dfy2];

 if(l>0)
 j=1;
 addx=[];
 addy=[];
 addd1=[];
 addd2=[];
 for i=1:m(1)
     if(Is_in_bound(lower_bound, upper_bound, x(i,:)))
        for k=1:narg
             addx(j,k)=x(i,k);
         end;
             addy(j)=y(i);
             addd1(j)=dfy1(i);
             addd2(j)=dfy2(i);
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
 y=[yy';addy';dfy1';addd1';dfy2';addd2'];
 else
     y=[yy';dfy1';dfy2'];
 end;


if (strcmp(str,'lin')==1)
    Xf=Polinom(1,x);
    X=Xf;
    for i=1:m(2)
        dX(:,:,i)=derivativePolinom(1,x,i);
        X=[X;dX(:,:,i)];
    end;
elseif (strcmp(str,'quad')==1)
    Xf=Polinom(2,x);
    X=Xf;
    for i=1:m(2)
        dX(:,:,i)=derivativePolinom(2,x,i);
        X=[X;dX(:,:,i)];
    end;
elseif (strcmp(str,'cub')==1)
    Xf=Polinom(3,x);
    X=Xf;
    for i=1:m(2)
        dX(:,:,i)=derivativePolinom(4,x,i);
        X=[X;dX(:,:,i)];
    end;
elseif (strcmp(str,'quar')==1)
    Xf=Polinom(4,x);
    X=Xf;
    dX=[];
    for i=1:m(2)
        dX(:,:,i)=derivativePolinom(4,x,i);
        X=[X;dX(:,:,i)];
    end;
end;


Xt=X';
A=Xt*X;
b=Xt*y;
a=A\b;


if (strcmp(str,'lin')==1)
    phio=@(x)sum(a'.*Polinom(1,x));
elseif (strcmp(str,'quad')==1)
    phio=@(x)sum(a'.*Polinom(2,x));
elseif (strcmp(str,'cub')==1)
    phio=@(x)sum(a'.*Polinom(3,x));
elseif (strcmp(str,'quar')==1)
    phio=@(x)sum(a'.*Polinom(4,x));
end;


    for i=1:length(r)
if(abs(upper_bound(i)-lower_bound(i))<eps)||(abs(coeff(i))<eps)
    break;
end;
if((upper_bound(i)<lower_bound(i)))
    swap(upper_bound(i),lower_bound(i));
end;
    end;

    
[r,out,opt]=fmincon(phio,rn,[],[],[],[],lower_bound,upper_bound)


l=l+1;
end;


disp('');
disp('');
disp('------------------------------RESULTS---------------------------------------------------------------');
disp('');
disp('');
xmin=r
fmin=f(xmin)
count=fInvokeCount

end

