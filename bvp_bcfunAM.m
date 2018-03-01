function resmu=bvp_bcfunAM(mu0,mue)
%function res=bcfun(y0,ye,params)
global params
bc_coef=params.bc;

switch params.type_bc(1)
    case 1
        resmu(1,1)=mu0(1);
    case 2
        pr=params.r;
        if strcmp(params.opt_u,'r')
            pr=Ut_fun(params.x0_xe(1));
        end
        resmu(1,1)=mu0(2)+(pr+bc_coef(1,2)/bc_coef(1,1))*mu0(1);  
end

switch params.type_bc(2)
    case 1
        resmu(2,1)=mue(1);
    case 2
        pr=params.r;
        if strcmp(params.opt_u,'r')
            pr=Ut_fun(params.x0_xe(2));
        end
        resmu(2,1)=mue(2)+(pr+bc_coef(2,2)/bc_coef(2,1))*mue(1);  
end

end