function yg=yg_fun(x)
% calculation yg (desirable value of y)

global params
nx=length(x);
yg=ones(1,nx);
% % yg=zeros(1,nx);
% % half_x0xe=0.5*(params.x0_xe(1)+params.x0_xe(2));
% % for i=1:nx
% %     if x(i)<=half_x0xe
% %         yg(i)=x(i);
% %     else
% %          yg(i)=params.x0_xe(2)-x(i);
% %     end
% % end