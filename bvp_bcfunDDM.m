function resz=bvp_bcfunDDM(z0,ze)
%function res=bcfun(y0,ye,params)
global params
bc_coef=params.bc;
resz=[bc_coef(1,1)*z0(2)+bc_coef(1,2)*z0(1)
         bc_coef(2,1)*ze(2)+bc_coef(2,2)*ze(1)];
end