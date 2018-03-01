function dmudx=bvp_odefunAM(x,mu)
%function dydx=odefun(x,y,params)
global params
global soly ib
%x
yg=yg_fun(x);
yx=deval(soly,x);
h_yg=2*params.gama_y*(yx(1)-yg);
Ut=Ut_fun(x);

switch params.opt_u
  case 'f0'
      dmudx=[mu(2)
                 -params.r*mu(2)+ params.g1*mu(1)+ 3*params.g3*yx(1)^2*mu(1)-h_yg];        
   case 'fu'
      dmudx=[mu(2)
                 -params.r*mu(2)+ params.g1*mu(1)+ 3*params.g3*yx(1)^2*mu(1)-h_yg];
             
   case 'r'
      dmudx=[mu(2)
                 -Ut*mu(2)+ params.g1*mu(1)+ 3*params.g3*yx(1)^2*mu(1)-h_yg];        

    case 'g1'
      dmudx=[mu(2)
                 -params.r*mu(2)+ Ut*mu(1)+ 3*params.g3*yx(1)^2*mu(1)-h_yg];        
 
    case 'g3'
       dmudx=[mu(2)
                 -params.r*mu(2)+ params.g1*mu(1)+ 3*Ut*yx(1)^2*mu(1)-h_yg]; 
         
end %switch params.opt_u

end