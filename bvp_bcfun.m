function resy=bvp_bcfun(y0,ye)
%function res=bcfun(y0,ye,params)
global params
bc_coef=params.bc;
% resy=[bc_coef(1,1)*y0(2)+bc_coef(1,2)*y0(1)-bc_coef(1,3)*bc_coef(1,4)
%          bc_coef(2,1)*ye(2)+bc_coef(2,2)*ye(1)-bc_coef(2,3)*bc_coef(2,4)];
resy(1,1)=[bc_coef(1,1)*y0(2)+bc_coef(1,2)*y0(1)-bc_coef(1,3)*bc_coef(1,4)];
resy(2,1) =[ bc_coef(2,1)*ye(2)+bc_coef(2,2)*ye(1)-bc_coef(2,3)*bc_coef(2,4)];

end