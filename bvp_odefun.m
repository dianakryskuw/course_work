function dydx=bvp_odefun(x,y)
%function dydx=odefun(x,y,params)
global params

% lts=params.lts;
% nu=floor(x/lts)+1;
% nu=min(nu,params.n);
% switch params.typeU
%     case 'const' 
%        Ut=params.U(nu);
%     case 'linear'
%         x1=params.x0_xe(1)+(nu-1)*lts;
%         x2=x1+lts;
%         Ut=(x2-x)/lts*params.U(nu)+(x-x1)/lts*params.U(nu+1);
%  
%     otherwise
%         error('Unexpected opproximation type - bvp_odefun')
% end
Ut=Ut_fun(x);
switch params.opt_u
  case 'f0'
      dydx=[y(2)
                 params.r*y(2)+ params.g1*y(1)+ params.g3*y(1)^3-Ut-params.fu];        
   case 'fu'
      dydx=[y(2)
                 params.r*y(2)+ params.g1*y(1)+ params.g3*y(1)^3-Ut-params.f0];  
   case 'r'
     dydx=[y(2)
               Ut*y(2)+ params.g1*y(1)+ params.g3*y(1)^3-params.f0-params.fu];
    case 'g1'
     dydx=[y(2)
                params.r*y(2)+ Ut*y(1)+ params.g3*y(1)^3-params.f0-params.fu];
    case 'g3'
     dydx=[y(2)
              params.r*y(2)+ params.g1*y(1)+ Ut*y(1)^3-params.f0-params.fu];
end %swich

%     dydx=[y(2)
%               params.r*y(2)+ params.g1*y(1)+ params.g3*y(1)^3-Ut-params.f];
%      dydx=[y(2)
%               -params.f];
% dydx = [y(2)
% -2*x*y(2) - 5*y(1) + cos(3*x)];
end