function psi0 = OptBVP_psi0(soly)
% Calculation of psi0 function

global params
global psi_y psi_u 

yx=deval(soly,params.x,1);
%psi_y=trapz(soly.x,(soly.y(1,:)-yg).^2);
psi_y=trapz(params.x,(yx-params.yg).^2);

psi_u=0;

if params.gama_u>0
 lts=params.lts;
 switch params.typeU
    case 'const' 
        for i=1:params.n
            psi_u=psi_u+params.U(i)^2;
        end
        psi_u=psi_u*lts;
   case 'linear'
        x1=params.x0_xe(1);
        for i=1:params.n
            x2=x1+lts;
            u1=params.U(i);
            u2=params.U(i+1);
            ai=x2*u1-x1*u2;
            bi=u2-u1;
            psi_u=psi_u+(ai^2*x2+ai*bi*x2^2+bi^2*x2^3/3);
            psi_u=psi_u-(ai^2*x1+ai*bi*x1^2+bi^2*x1^3/3);
            x1=x2;
        end
      psi_u=psi_u/lts^2;
     otherwise
        error('Unexpected opproximation type /odefun')
 end   %switch params.typeU
end   %if params.gama>0

psi0 =params.gama_y*psi_y+params.gama_u*psi_u;

end

