function dzdx=bvp_odefunDDM(x,z)
%function dydx=odefun(x,y,params)
global params
global soly ib

du_db=dudb_fun(x,ib);
% % lts=params.lts;
% % 
% % nu=floor(x/lts)+1;
% % nu=min(nu,params.n);
% % du_db=0;
% % switch params.typeU
% %     case 'const' 
% %         if nu==ib
% %             du_db=1;
% %         end %if nu=ib
% %     case 'linear'
% %         x1=params.x0_xe(1)+(nu-1)*lts;
% %         x2=x1+lts;
% %         %
% %        % Ut=(x2-x)/lts*params.U(nu)+(x-x1)/lts*params.U(nu+1);
% %         
% %         switch ib
% %             case 1
% %                 if nu==1
% %                    du_db=(x2-x)/lts ;
% %                 end % if nu==1 
% %                  
% %               case params.n+1
% %                 if nu==params.n
% %                    du_db=(x-x1)/lts ;
% %                 end % if nu==n
% %                 
% %               otherwise
% %                  if nu==ib-1
% %                    du_db=(x-x1)/lts;
% %                 end % if nu==ib-1
% %                  if nu==ib
% %                    du_db=(x2-x)/lts;
% %                 end % if nu                  
% %         end %switch ib
% %   
% %     otherwise
% %         error('Unexpected opproximation type  - bvp_odefunDDM')
% % end %switch params.typeU

yx=deval(soly,x);
switch params.opt_u
  case 'f0'
      dzdx=[z(2)
                 params.r*z(2)+ params.g1*z(1)+ 3*params.g3*yx(1)^2*z(1)-du_db];        
   case 'fu'
       dzdx=[z(2)
                 params.r*z(2)+ params.g1*z(1)+ 3*params.g3*yx(1)^2*z(1)-du_db];        
 
   case 'r'
      dzdx=[z(2)
                 params.r*z(2)+du_db* yx(2)+params.g1*z(1)+ 3*params.g3*yx(1)^2*z(1)];        

    case 'g1'
      dzdx=[z(2)
                 params.r*z(2)+du_db* yx(1)+params.g1*z(1)+ 3*params.g3*yx(1)^2*z(1)];        
 
    case 'g3'
      dzdx=[z(2)
                 params.r*z(2)+params.g1*z(1)+ 3*params.g3*yx(1)^2*z(1)+du_db* yx(1)^3];        
        
end %switch params.opt_u

end