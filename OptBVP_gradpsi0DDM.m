function grad_psi0_ib= OptBVP_gradpsi0DDM(solz)
% Calculation of partial derivative of psi0 
%with respecr to ib component of optimization vector

global params
global soly ib

%x=linspace(params.x0_xe(1),params.x0_xe(2),params.n_gradpsi0);
yx=deval(soly,params.x,1);
zx=deval(solz,params.x,1);

grad_psi0_y=trapz(params.x,(yx-params.yg).*zx);
%grad_psi_u=0;
grad_psi_u=trapz(params.x,Ut_fun(params.x).*dudb_fun(params.x,ib));
% if params.gama>0
%  lts=params.lts;
%  switch params.typeU
%     case 'const' 
%         for i=1:params.n
%             psi_u=psi_u+params.U(i)^2;
%         end
%         psi_u=psi_u*lts;
%    case 'linear'
%         x1=params.x0_xe(1);
%         for i=1:params.n
%             x2=x1+lts;
%             u1=params.U(i);
%             u2=params.U(i+1);
%             ai=x2*u1-x1*u2;
%             bi=u2-u1;
%             psi_u=psi_u+(ai^2*x2+ai*bi*x2^2+bi^2*x2^3/3);
%             psi_u=psi_u-(ai^2*x1+ai*bi*x1^2+bi^2*x1^3/3);
%             x1=x2;
%         end
%       psi_u=psi_u/lts^2;
%      otherwise
%         error('Unexpected opproximation type /odefun')
%  end   %switch params.typeU
% end   %if params.gama>0

grad_psi0_ib =2*params.gama_y*grad_psi0_y+2*params.gama_u*grad_psi_u;
%psi0=0.5*trapz(sol.x,(sol.y(1,:)-params.yg).^2)   

end