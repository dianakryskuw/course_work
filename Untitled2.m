function [  ] = Untitled2( str ,np,bound1,bound2,bound3,startPoint )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


global params
global psi_y psi_u 
global k_func k_grad
k_func=0; k_grad=0;

%interval
x0=0;
xe=1;
params.x0_xe=[x0 xe];

%coefficients of equation
        f0=1;
        fu=3;
        p=1;
        r=5;
        g1=6;
        g3=7;

% % f0=0;
% % fu=0;
% % p=0;
% % r=0;
% % g1=1;
% % g3=0;
params.f0=f0;
params.fu=fu;
params.p=p;
params.r=r;
params.g1=g1;
params.g3=g3;

%boundary conditions
%bc at  point x0: a0*dy/dx+b0*y=d0*yg_x0
%     at  point xe: ae*dy/dx+be*y=de*yg_x0
%  bc_coef=[a0 b0 d0 yg_x0;
%                 ae be de yg_xe];
%  a0,ae=1 for Neumann boundary condition
%  a0,ae=0 for Dirichlet boundary condition

%bc_coef=[0 1 1 0;0 1 1  0]; % Dirichlet boundary condition
bc_coef=[1 0 1 1;1 0 1  -1]; % Neumann boundary condition
params.bc=bc_coef;

% Dirichlet boundary condition on x0 and xe
%params.type_bc=[1 1]; % Dirichlet boundary condition
% % Neumann boundary condition on x0 and xe
params.type_bc=[2 2]; % Neumann boundary condition

n=3;%number of interval !!!
params.n=n;
lx=xe-x0;
%length of interval !!! All interval are equal
params.lts=lx/n;

% beta=0.5;%???
% params.beta=beta;

        %Optimization parameters
       %params.opt_u='f0'
       %params.opt_u='fu'
       params.opt_u='r'
     %params.opt_u='g1'
      %params.opt_u='g3'
switch params.opt_u
  case 'f0'
      u=f0;
   case 'fu'
      u=fu;  
   case 'r'
      u=r;
    case 'g1'
      u=g1;
    case 'g3'
      u=g3;
end %swich

%Type of optimization function
typeU='const'
%typeU='linear';
params.typeU=typeU;

%Initial values of optimization parameters
switch params.typeU
    case 'const' 
       nu=n;
       U0=u*ones(1,nu)
      case 'linear'
       nu=n+1;
       U0=u*ones(1,nu)
    otherwise
        error('Unexpected opproximation type')
end

% U0=[0.0110    0.6025    4.9985    0.3861    0.0393]

params.U=U0;

% Lower and upper bounds on the optimization parameters
% % un=U0-alf*abs(U0);
% % uv=U0+alf*abs(U0);

un=-10*ones(1,nu)
uv=10*ones(1,nu)
params.un=un;
params.uv=uv;

% Number of point for calculation  of optimization criteria 
% of integral type
n_gradpsi0=50; 
params.n_gradpsi0=n_gradpsi0;

% Multipliers for calculation  of optimization criteria 
gama_y=1;
gama_u=1e-2;
%gama_u=0;
params.gama_y=gama_y;
params.gama_u=gama_u;

% Desirable values of BVP solution
yg=ones(1,params.n_gradpsi0);
       x=linspace(params.x0_xe(1),params.x0_xe(2),n_gradpsi0);
        params.x=linspace(params.x0_xe(1),params.x0_xe(2),params.n_gradpsi0);
% %         half_x0xe=0.5*(params.x0_xe(1)+params.x0_xe(2));
% %         for i=1:n_gradpsi0
% %             if x(i)<=half_x0xe
% %                 yg(i)=x(i);
% %             else
% %                  yg(i)=params.x0_xe(2)-x(i);
% %             end
% %         end
params.yg=yg;

% plot U0
figure 
Ut=Ut_fun(x);
hold on
plot(x,Ut,'b')
title('') 



% figure
% plot(params.x,params.yg)
% title('yg')
% p1=0.25;
% p2=0.75;
% x1=x0+p1*lx;
% x2=x0+p2*lx;


% params.y1d=y1d;
% params.y2d=y2d;
% params.y1_ad=y1_ad;
% params.y2_ad=y2_ad;


% method of derivative calculation
% 'FDM' - finite difference method 
% 'DDM' - direct difference method 
%'AM' - finite difference method 
params.deriv='DDM';

% Number of points in which we guess initial values 
% of BVP problem solution
ninit=30;
params.ninit=ninit;

params.plot='Yes';

figure
plot(params.x,params.yg,'b')

params.sol_color='m';
params.plot='Yes';
OptBVP_FuncDeriv([5,5,5]);
params.plot='No';

params.sol_color='r';
params.plot='No';

tic
[opx,opf,ope,opout]= fmincon(@OptBVP_FuncDeriv,[6,6,6],[],[],[],[],[-10,-10,-10],[10,10,10])
toc
figure(1)
 %  x=linspace(params.x0_xe(1),params.x0_xe(2),n_gradpsi0);
Ut=Ut_fun(x);
plot(x,Ut,'r')

params.plot='Yes';
figure(2)
OptBVP_FuncDeriv(opx);
params.plot='No';
x1l=bound1(1);
x1u=bound1(2);
x2l=bound2(1);
x2u=bound2(2);
 x3l=bound3(1);
 x3u=bound3(2);
%  x4l=bound4(1);
%  x4u=bound4(2);
Z1=x1l:0.1:x1u;
Z2=x2l:0.1:x2u;
[z1,z2]=meshgrid(Z1,Z2);

%f=str2func(fun);


r=startPoint;
%rn=[r(1)+100,r(2)+100];
rn=[r(1)+100,r(2)+100,r(3)+100];
%rn=[r(1)+100,r(2)+100,r(3)+100,r(4)+100];
eps=0.001;
zz=[z1,z2];
%mesh(z1,z2,f(z));
hold on
x1=[];
addx=[];
addy=[];
y=[];
n=np;
l=0;
k1=0.2*(x1u-x1l);
k2=0.2*(x2u-x2l);
k3=0.2*(x3u-x3l);
%k4=0.2*(x4u-x4l);
fInvokeCount=0;
tic

while ((l==0)||(norm(rn-r)>1)||(abs(phio(r)-OptBVP_FuncDeriv(r))>eps))
    if(l>0)
    k1=abs(rn(1)-r(1))*1.5;
    k2=abs(rn(2)-r(2))*1.5;
    k3=abs(rn(3)-r(3))*1.5;
   % k4=abs(rn(4)-r(4))*1.5;
    end;
    rn=r;
a1=max(r(1)-k1,x1l);
a2=min(r(1)+k1,x1u);
b1=max(r(2)-k2,x2l);
b2=min(r(2)+k2,x2u);
c1=max(r(3)-k3,x3l);
c2=min(r(3)+k3,x3u);

%d1=max(r(4)-k4,x4l);
%d2=min(r(4)+k4,x4u);


%plot3(r(1),r(2),f(r(1),r(2)),'*r');
x=transf(a1,a2,b1,b2,c1,c2,n);
%x=transf2(a1,a2,b1,b2,n);
%x=transf4(a1,a2,b1,b2,c1,c2,d1,d2,n);

j=1;
addx=[];
addy=[];
for i=1:length(x1)
    if((x1(i)>a1)&&(x1(i)<a2)&&(x2(i)>b1)&&(x2(i)<b2)&&(x3(i)>c1)&&(x3(i)<c2))
        addx(j,1)=x1(i);
        addx(j,2)=x2(i);
        addx(j,3)=x3(i);
%        addx(j,4)=x4(i);
        addy(j)=y(i);
        j=j+1;
    end;
end;
x1d=x(:,1);
x2d=x(:,2);
x3d=x(:,3);
%x4d=x(:,4);

for i=1:length(x1d)
        %myx=[x(i,1),x(i,2),x(i,3),x(i,4)];
        myx=[x(i,1),x(i,2),x(i,3)];
       % myx=[x(i,1),x(i,2)];
    yy(i)=OptBVP_FuncDeriv(myx);
    fInvokeCount=fInvokeCount+1;
    %yy(i)=f(x1d(i),x2d(i));
    %yy(i)=100*(x2d(i)-x1d(i)^2)^2+(1-x1d(i))^2;
    %yy(i)=x1(i).^4+x2(i).^4-10;
%yy(i)=2.*x1(i).^2-4.*x1(i)+x2(i).^2-8.*x2(i)+3;
end;


if(~isempty(addx))
x1=[x1d;addx(:,1)];
x2=[x2d;addx(:,2)];
x3=[x3d;addx(:,3)];
%x4=[x4d;addx(:,4)];
else
x1=x1d;
x2=x2d;
x3=x3d;
%x4=x4d;
end;
if(~isempty(addy))
y=[yy';addy'];
else
    y=yy';
end;


if (strcmp(str,'lin')==1)
    X=[ x1.^0, x1, x2,x3];
elseif (strcmp(str,'quad')==1)
    
%    X=[ x1.^0, x1, x2,x3,x4, x1.*x2, x1.*x3,x1.*x4,x2.*x3,x2.*x4,x3.*x4,x1.^2, x2.^2,x3.^2,x4.^2];
     X=[ x1.^0,  x1, x2,x3, x1.*x2, x1.*x3,x2.*x3,x1.^2, x2.^2,x3.^2];
    %X=[ x1.^0,  x1, x2, x1.*x2, x1.^2, x2.^2];
elseif (strcmp(str,'cub')==1)
    X=[ x1.^0,  x1, x2,x3, x1.*x2, x1.*x3,x2.*x3,x1.^2, x2.^2,x3.^2, x1.*x2.^2,x1.*x3.^2,x2.*x1.^2,x2.*x3.^2,x3.*x1.^2,x3.*x2.^2,x1.^3,x2.^3,x3.^3];
elseif (strcmp(str,'quar')==1)
    X=[ x1.^0, x1, x2, x1.*x2, x1.^2, x2.^2, x1.*x2.^2, x1.*x2.^2, x1.^3, x2.^3, x1.*x2.^3, x2.*x1.^3, x1.^2.*x2.^2, x1.^4, x2.^4];
end;


Xt=X';
A=Xt*X;
b=Xt*y;
a=A\b;


if (strcmp(str,'lin')==1)
    phi=@(x1,x2,x3)a(1)+a(2).*x1+a(3).*x2+a(4).*x3;
    phio=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(3);
elseif (strcmp(str,'quad')==1)
    
    
   % phio=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(3)+a(5).*x(4)+a(6).*x(1).*x(2)+a(7).*x(1).*x(3)+a(8).*x(1).*x(4)+a(9).*x(2).*x(3)+a(10).*x(2).*x(4)+a(11).*x(3).*x(4)+a(12).*x(1).^2+a(13).*x(2).^2+a(14).*x(3).^2+a(15).*x(4).^2;
    
  %  phi=@(x1,x2,x3)a(1)+a(2).*x1+a(3).*x2+a(4).*x3+a(5).*x1.*x2+a(6).*x1.*x3+a(7).*x2.*x3+a(8).*x1.^2+a(9).*x2.^2+a(10).*x3.^2;
    phio=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(3)+a(5).*x(1).*x(2)+a(6).*x(1).*x(3)+a(7).*x(2).*x(3)+a(8).*x(1).^2+a(9).*x(2).^2+a(10).*x(3).^2;
    %phi=@(x1,x2)a(1)+a(2).*x1+a(3).*x2+a(4).*x1.*x2+a(5).*x1.^2+a(6).*x2.^2;
    %phio=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(1).*x(2)+a(5).*x(1).^2+a(6).*x(2).^2;
elseif (strcmp(str,'cub')==1)
    phi=@(x1,x2,x3)a(1)+a(2).*x1+a(3).*x2+a(4).*x3+a(5).*x1.*x2+a(6).*x1.*x3+a(7).*x2.*x3+a(8).*x1.^2+a(9).*x2.^2+a(10).*x3.^2+a(11).*x1.*x2.^2+a(12).*x1.*x3.^2+a(13).*x2.*x1.^2+a(14).*x2.*x3.^2+a(15).*x3.*x1.^2+a(16).*x3.*x2.^2+a(17).*x1.^3+a(18).*x2.^3+a(19).*x3.^3;
   phio=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(3)+a(5).*x(1).*x(2)+a(6).*x(1).*x(3)+a(7).*x(2).*x(3)+a(8).*x(1).^2+a(9).*x(2).^2+a(10).*x(3).^2+a(11).*x(1).*x(2).^2+a(12).*x(1).*x(3).^2+a(13).*x(2).*x(1).^2+a(14).*x(2).*x(3).^2+a(15).*x(3).*x(1).^2+a(16).*x(3).*x(2).^2+a(17).*x(1).^3+a(18).*x(2).^3+a(19).*x(3).^3;
elseif (strcmp(str,'quar')==1)
    phi=@(x1,x2)a(1)+a(2).*x1+a(3).*x2+a(4).*x1.*x2+a(5).*x1.^2+a(6).*x2.^2+a(7).*x2.*x1.^2+a(8).*x1.*x2.^2+a(9).*x1.^3+a(10).*x2.^3+a(11).*x1.*x2.^3+a(12).*x2.*x1.^3+a(13).*x1.^2.*x2.^2+a(14).*x1.^4+a(15).*x2.^4;
    phio=@(x)a(1)+a(2).*x(1)+a(3).*x(2)+a(4).*x(1).*x(2)+a(5).*x(1).^2+a(6).*x(2).^2+a(7).*x(2).*x(1).^2+a(8).*x(1).*x(2).^2+a(9).*x(1).^3+a(10).*x(2).^3+a(11).*x(1).*x(2).^3+a(12).*x(2).*x(1).^3+a(13).*x(1).^2.*x(2).^2+a(14).*x(1).^4+a(15).*x(2).^4;
end;


if((abs(a1-a2)<eps)||(abs(b1-b2)<eps) || (abs(k1)<eps) || (abs(k2)<eps))
    break;
end;
if(a1>a2)
ta=a1;
a1=a2;
a2=ta;
% k=k+0.02;
end;

if(b1>b2)
tb=b1;
b1=b2;
b2=tb;
% k=k+0.02;
 end;

%[r,out,opt]=fmincon(phio,rn,[],[],[],[],[a1,b1],[a2,b2])
 [r,out,opt]=fmincon(phio,rn,[],[],[],[],[a1,b1,c1],[a2,b2,c2])
 % [r,out,opt]=fmincon(phio,rn,[],[],[],[],[a1,b1,c1,d1],[a2,b2,c2,d2])


if(norm(rn-r)<eps)
    break;
end;


%f(r(1),r(2))
%if(abs(phio(r)-f(r(1),r(2)))<1)
%      k1=abs(rn(1)-r(1))
%      k2=abs(rn(2)-r(2))
%end;
l=l+1;
ff(l)=OptBVP_FuncDeriv(r);


end;
toc
disp('');
disp('');
disp('------------------------------RESULTS---------------------------------------------------------------');
disp('');
disp('');
xmin=r
fmin=OptBVP_FuncDeriv(xmin)


fInvokeCount
ff
%plot3(r(1),r(2),f(r(1),r(2)),'ok');

%mesh(z1,z2,phi(z1,z2));
 k_func
 k_grad
params.plot='Yes'; 
hold on
params.sol_color='g';
[f,g]=OptBVP_FuncDeriv(xmin)

grid on
disp('psi_y psi_u')
psi_y
psi_u

figure(1)
       x=linspace(params.x0_xe(1),params.x0_xe(2),n_gradpsi0);
Ut=Ut_fun(x);
plot(x,Ut,'g')
grid on

return