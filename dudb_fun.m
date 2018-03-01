function du_db=dudb_fun(x,ib)
global params
%global soly ib

lts=params.lts;
nx=length(x);
du_db=zeros(1,nx);

for i=1:nx
 nu=floor(x(i)/lts)+1;
 nu=min(nu,params.n);
 switch params.typeU
    case 'const' 
        if nu==ib
            du_db(i)=1;
        end %if nu=ib
    case 'linear'
        x1=params.x0_xe(1)+(nu-1)*lts;
        x2=x1+lts;
        %
        %Ut=(x2-x)/lts*params.U(nu)+(x-x1)/lts*params.U(nu+1);
        
        switch ib
            case 1
                if nu==1
                   du_db(i)=(x2-x(i))/lts ;
                end % if nu==1 
                 
              case params.n+1
                if nu==params.n
                   du_db(i)=(x(i)-x1)/lts ;
                end % if nu==n
                
              otherwise
                 if nu==ib-1
                   du_db(i)=(x(i)-x1)/lts;
                end % if nu==ib-1
                 if nu==ib
                   du_db(i)=(x2-x(i))/lts;
                end % if nu                  
        end %switch ib
  
    otherwise
        error('Unexpected opproximation type  - bvp_odefunDDM')
 end %switch params.typeU
end %for i=1:nx
end
