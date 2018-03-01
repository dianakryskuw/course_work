function [c,ceq] = optconstr1(u,params)
c = [];

params.U=u;
options = odeset('RelTol',1e-5,'AbsTol',1e-5);
sol = ode15s(@LV_2r4p_F, params.t0_te, params.Y0, options,params); 

%psi0= LV_2r4p_psi0(sol,params);

switch params.n_Opt_kr
    case 1
     %variant 5a n_Opt_par 1-3 
     a=sol.y(2,:)-params.y2_ad;
     a=a+abs(a);
     ceq =trapz(sol.x,a.^2); 
  
    case 2
     %variant 5b n_Opt_par 4-7 
     a=sol.y(1,:)-params.y1_ad;
     a=a+abs(a);
     ceq =trapz(sol.x,a.^2); 
     
end


