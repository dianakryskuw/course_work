function [psi0,grad_psi0]= OptBVP_FuncDeriv(u)
% function [psi0,grad_psi0]= OptFuncDeriv(u,params)

global params
global soly ib
global k_func k_grad

k_func=k_func+1;

params.U=u;

xinit=linspace(params.x0_xe(1),params.x0_xe(2),params.ninit);
yinit=[10 0];
solinit=bvpinit(xinit,yinit);
% solinit.x
% solinit.y

soly=bvp5c(@bvp_odefun,@bvp_bcfun,solinit);
%disp('OptBVP_FuncDeriv');
psi0= OptBVP_psi0(soly);

%plot solution
if strcmp(params.plot,'Yes')
   % figure
   hold on
    yx=deval(soly,params.x,1);
    plot(params.x, yx,params.sol_color)
    %plot(soly.x, soly.y(1,:),'r')
%     xlabel('x') ; ylabel('y')
%     hold on
  %  plot(params.x,params.yg,params.sol_color)
end

 if nargout>1
     switch params.typeU
      case 'const' 
       nb=params.n;
      case 'linear'
       nb=params.n+1;
     otherwise
        error('Unexpected opproximation type - /OptBVP_FuncDeriv/')
     end %switch params.typeU

  grad_psi0=zeros(1,nb);
  switch params.deriv
    case 'FDM'
      %finite difference method (FDM)
      disp('finite difference method');

  
     disp('____  end finite difference method _____');
    
    case 'DDM'
        k_grad=k_grad+1;
      % direct difference method (DDM)
      disp('direct difference method');
          %xb=linspace(params.x0_xe(1),params.x0_xe(2),20)
        for ib=1:nb
            xinit=linspace(params.x0_xe(1),params.x0_xe(2),params.ninit);
            zinit=[10 0];
            solinit=bvpinit(xinit,zinit);
            solz=bvp5c(@bvp_odefunDDM,@bvp_bcfunDDM,solinit);
            
            grad_psi0(ib)= OptBVP_gradpsi0DDM(solz);
           
        end % for ib=1:nb
      disp('___ end direct difference method ___');
      
    case 'AM'
          k_grad=k_grad+1;
        disp('adjoint method (AM)')
           disp('adjoint method');
            xinit=linspace(params.x0_xe(1),params.x0_xe(2),params.ninit);
            zinit=[10 0];
            solinit=bvpinit(xinit,zinit);
            solmu=bvp5c(@bvp_odefunAM,@bvp_bcfunAM,solinit);
            
            grad_psi0= OptBVP_gradpsi0AM(solmu,nb);      
      
          disp('___ end adjoint method ___');      
    end %switch
 end %if nargout>1
end %function LVpar_OptFuncDeriv4



