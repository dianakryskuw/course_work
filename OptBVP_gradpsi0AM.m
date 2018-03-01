function grad_psi0= OptBVP_gradpsi0AM(solmu,nb)
% Calculation of partial derivative of psi0 
%with respecr to ib component of optimization vector

global params
global soly ib

%x=linspace(params.x0_xe(1),params.x0_xe(2),params.n_gradpsi0);
mux=deval(solmu,params.x,1);

grad_psi0=zeros(1,nb);
gu=2*params.gama_u*Ut_fun(params.x);
switch params.opt_u
    case 'f0'
        for i=1:nb
          %gu=2*params.gama_u*Ut_fun(params.x);
          grad_psi0(i) =trapz(params.x,(gu+mux).*dudb_fun(params.x,i));  
        end  

     case 'fu'
        for i=1:nb
          %gu=2*params.gama_u*Ut_fun(params.x);
          grad_psi0(i) =trapz(params.x,(gu+mux).*dudb_fun(params.x,i));  
        end  
             
    case 'r'
        yx=deval(soly,params.x,2);
        for i=1:nb
          gu=2*params.gama_u*Ut_fun(params.x);
          grad_psi0(i) =trapz(params.x,(gu-mux.*yx).*dudb_fun(params.x,i));  
        end  
        
    case 'g1'
        yx=deval(soly,params.x,1);
        for i=1:nb
          %gu=2*params.gama_u*Ut_fun(params.x);
          grad_psi0(i) =trapz(params.x,(gu-mux.*yx).*dudb_fun(params.x,i));  
        end  
 
    case 'g3'
        yx=deval(soly,params.x,1).^3;
        for i=1:nb
         % gu=2*params.gama_u*Ut_fun(params.x);
          grad_psi0(i) =trapz(params.x,(gu-mux.*yx).*dudb_fun(params.x,i));  
        end   
         
end %switch params.opt_u

end