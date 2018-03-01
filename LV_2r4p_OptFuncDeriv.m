function [psi0,grad_psi0]= LV_2r4p_OptFuncDeriv(u,params)
%function [psi0,grad_psi]= LV_2r4p_OptFuncDeriv(u,params)
% LVpar_OptFuncDeriv4 calculates fuction f, which is criteria 
% for identification parameters problem and derivative g=df/du
%  for the following ODE Volterra model 
%    y(1)'=p(1,1)*y1+p(1,2)*y(1)*y(2)+...+p(1,kd)*...
%     . . . . . . . .. . . . . . . . . . . . . 
%    y(krv)'=p(kr,1)*y(krv)+p(kr,2)*y(1)*y(kr)+...+p(kr,kd)*...
%
% Input parameters:
%   u -values of unknown parametes p(i,j);
%   params.kr - number of equations in Volterra model 
%   params.kd - number of term in each equation
%
%              |p(1,1)   ...   p(1,kd)  |   
%   params.p=  |.  .  .  .  .  .  .  .   | - model parameters
%              |p(kr,1)  ...  p(krv,kd)|   
%
%               |indu(1,1)  ... indu(1,kd) |   indicator of 
%   params.indu=|.  .  .  .  .  .  .   .   | - identification
%               |indu(kr,1) ... indu(kr,kd)|   model parameters,
%
%      indu(i,j)=1 - parameters p(i,j) is subject to uidentification
%      indu(i,j)=0 - parameters p(i,j) is known
%
%   params.t0_te=[t0 te] - interval of integration of ODE
%   params.Y0 - initial conditions
%   params.plot - const for plotting solution and diravative
%                ('Yes' - solution and diravative are plotting 
%                 ('No' - solution and diravative are not plotting)
%   params.deriv - indicator of metod for calcution of derivative
%                  ('FDM' - finite difference method;
%                   'DDM' -direct difference method )
%   params.ph - matrix of increment of p parameters when FDM is used for 
%         calculation of derivative matrix g
%
%   params.tz - coordinates of points in which solution y1, y2 is known
%   params.yz=[y1z; y2z] - known solution y1, y2 in points tz
%
% Output parameters:psi0,grad_psi0


params.U=u;
options = odeset('RelTol',1e-5,'AbsTol',1e-5);
sol = ode15s(@LV_2r4p_F, params.t0_te, params.Y0, options,params); 

psi0= LV_2r4p_psi0(sol,params);

%plot solution
if strcmp(params.plot,'Yes')
   figure
   plot(sol.y(1, :), sol.y(2,:))
   title('solution y1(t),y2(t)')
   xlabel('y1')
   ylabel('y2')
   grid 'on'
   figure
   plot(sol.x,sol.y(1,:))
   title('solution t,y1(t)')
   xlabel('t')
   ylabel('y1')
   grid 'on'
   figure
   plot(sol.x,sol.y(2,:))
   title('solution t, y2(t)')
   xlabel('t')
   ylabel('y2')
   grid 'on'
   nu=params.n+1
   tu=zeros(1,nu);
   for i=1:nu
       tu(i)=params.t0_te(1)+(i-1)*params.lts;
   end
   tu
  figure
  plot(tu,u,tu,params.un,tu,params.uv)
  xlabel('t')
  ylabel('u')
   grid 'on'
   u
   disp('params.p')
   disp(params.p)
   psi0
end %if strcmp(params.plot,'Yes')

if nargout==1
 return
end
n=params.n;
g=zeros(1,n);
grad_psi0=g;
switch params.deriv
  case 'FDM'
   %finite difference method (FDM)
   disp('finite difference method');

  
    disp('____  end finite difference method _____');
    
    
  case 'DDM'
    %direct difference method (DDM)
    disp('direct difference method');
    options1 = odeset('RelTol',1e-5,'AbsTol',1e-5);
    Y0D=[0 0];
    %m2=length(params.tz);
    
 
    disp('___ end direct difference method ___');
end %switch
end %function LVpar_OptFuncDeriv4
