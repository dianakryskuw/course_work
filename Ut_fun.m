function Ut=Ut_fun(x)
%Calculation of optimization parameters values
%in points x

global params

lts=params.lts;
nx=length(x);
Ut=zeros(1,nx);

for i=1:nx
 nu=floor(x(i)/lts)+1;
 nu=min(nu,params.n);
 switch params.typeU
    case 'const' 
       Ut(i)=params.U(nu);
    case 'linear'
        x1=params.x0_xe(1)+(nu-1)*lts;
        x2=x1+lts;
        Ut(i)=(x2-x(i))/lts*params.U(nu)+(x(i)-x1)/lts*params.U(nu+1);
     otherwise
        error('Unexpected opproximation type - bvp_odefun')
 end % switch params.typeU
end % for i=1:nx
end
