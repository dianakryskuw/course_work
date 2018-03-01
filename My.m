function [ f ] = My(testPoints)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% clc
% clear
% clear global 
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
       U0=u*ones(1,nu);
      case 'linear'
       nu=n+1;
       U0=u*ones(1,nu);
    otherwise
        error('Unexpected opproximation type')
end

U0=testPoints
% U0=[0.0110    0.6025    4.9985    0.3861    0.0393]

params.U=U0

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
title('optimization function - initial (b), optimal (r)') 



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

params.sol_color='r';
params.plot='No';
            f=OptBVP_FuncDeriv(testPoints)
             disp('psi_y psi_u')
             psi_y
             psi_u
            return

tic
[f,g]=OptBVP_FuncDeriv(U0)
%title('yg (blue),  initial solution y (red)')
toc
grid on
disp('psi_y psi_u')
psi_y
psi_u


params.deriv='AM';

% Number of points in which we guess initial values 
% of BVP problem solution
ninit=30;
params.ninit=ninit;

params.plot='Yes';

%f=OptBVP_FuncDeriv(U0)
tic
[f,g]=OptBVP_FuncDeriv(U0)
toc
%title('yg (blue),  init so y (red)')
%grid on
disp('psi_y psi_u')
psi_y
psi_u

%[c,ceq] = optconstr1(U0,params)
%options1=optimset('MaxFunEvals',1100);
params.deriv='AM';
params.plot='No';
x0 = U0
k_func=0; k_grad=0;

options = optimset('GradObj','on','LargeScale','on');
%options = optimset('GradObj','off','LargeScale','on','MaxFunEvals',0,'MaxIter',0);
tic
[u_opt,f_opt,exitflag,output,lambda,grad,hessian] = fmincon(@OptBVP_FuncDeriv,x0,...
         [],[],[],[],un,uv,[],options)
 toc
  
 k_func
 k_grad
params.plot='Yes'; 
hold on
params.sol_color='g';
[f,g]=OptBVP_FuncDeriv(u_opt)

title('yg (blue),  init  y (red), opt y (green)')
grid on
disp('psi_y psi_u')
psi_y
psi_u

figure(1)
Ut=Ut_fun(x);
plot(x,Ut,'r')
grid on

% [c,ceq] = optconstr1(x,params)
% disp('psi_y psi_u')
% psi_y
% psi_u
% % params.plot='No';
% % x0 = U0
% % % [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@LV_2r4p_OptFuncDeriv,x0,...
% % %          [],[],[],[],un,uv,@optconstr1,options,params)
% % [x,fval,exitflag,output] = fmincon(@LV_2r4p_OptFuncDeriv,x0,...
% %          [],[],[],[],un,uv,@optconstr1,options,params)
% % 
% % params.plot='Yes';
% % f=LV_2r4p_OptFuncDeriv(x,params)
% % [c,ceq] = optconstr1(x,params)
% % %end
end

