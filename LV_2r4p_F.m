function F = LV_2r4p_F(t, y, params)
%calculation of ODE right-hand

n_Opt_par=params.n_Opt_par;
lts=params.lts;
nu=floor(t/lts)+1;
nu=min(nu,params.n);
switch params.typeU
    case 'const' 
       Ut=params.U(nu);
    case 'linear'
        t1=params.t0_te(1)+(nu-1)*lts;
        t2=t1+lts;
        Ut=(t2-t)/lts*params.U(nu)+(t-t1)/lts*params.U(nu+1);
 
    otherwise
        error('Unexpected opproximation type')
end
p=params.p;
p(n_Opt_par)=Ut;


%         dy1dt: 'p1*y1 + p2*y1^2 - p5*y1*y2*(exp(-p6*y1) - 1)'
%         dy2dt: 'p3*y2 + p4*y2^2 - p5*p7*y1*y2*(exp(-p6*y1) - 1)'

% F = [p(1)*y(1)+ p(2)*y(1)^2 - p(5)*y(1)*y(2)*(exp(-p(6)*y(1)) - 1)
%      p(3)*y(2) + p(4)*y(2)^2 - p(5)*p(7)*y(1)*y(2)*(exp(-p(6)*y(1)) - 1)];
 

 switch params.model
    case 1
%         dy1dt: 'p1*y1 + p2*y1^2 - p5*y1*y2*(exp(-p6*y1) - 1)'
%         dy2dt: 'p3*y2 + p4*y2^2 - p5*p7*y1*y2*(exp(-p6*y1) - 1)'
         A1=p(1)*y(1)+ p(2)*y(1)^2;
         B1=-p(5)*y(1)*y(2)*(exp(-p(6)*y(1)) - 1);
         A2=p(3)*y(2) + p(4)*y(2)^2;
         B2=p(7)*B1;
         F=[A1+B1;A2+B2];

%      F = [p(1)*y(1)+ p(2)*y(1)^2 - p(5)*y(1)*y(2)*(exp(-p(6)*y(1)) - 1)
%          p(3)*y(2) + p(4)*y(2)^2 - p(5)*p(7)*y(1)*y(2)*(exp(-p(6)*y(1)) - 1)];
    case 2
%         dy1dt: 'p1*y1 + (p4*y1^2*y2)/((p5*y1^2 + 1)*(p6 + p7*y2))'
%         dy2dt: 'p2*y2+p3*y2^2+(p4*p8*y1^2*y2^2)/((p5*y1^2+1)*(p6+p7*y2)*(p9+p10*y2))'
        A1=p(1)*y(1);
        B1=p(4)*y(1)^2*y(2)/((p(5)*y(1)^2 + 1)*(p(6) + p(7)*y(2)));
        A2=p(2)*y(2)+p(3)*y(2)^2;
        B2=p(8)*y(2)/(p(9)+p(10)*y(2))*B1;
        F=[A1+B1;A2+B2];
%         
%        F=[p(1)*y(1) + p(4)*y(1)^2*y(2)/((p(5)*y(1)^2 + 1)*(p(6) + p(7)*y(2)))
%         p(2)*y(2)+p(3)*y(2)^2+p(4)*p(8)*y(1)^2*y(2)^2/((p(5)*y(1)^2+1)*(p(6)+p(7)*y(2))*(p(9)+p(10)*y(2)))];

    case 3
%         dy1dt: 'p1*y1 - (p4*y1*y2*(exp(-p5*y1) - 1))/(p6 + p7*y2)'
%         dy2dt: 'p2*y2 + p3*y2^2 - (p4*p8*y1*y2*(exp(-p5*y1) - 1))/ (p6 + p7*y2)'
        A1=p(1)*y(1);
        B1=-(p(4)*y(1)*y(2)*(exp(-p(5)*y(1)) - 1))/(p(6) + p(7)*y(2));
        A2=p(2)*y(2) + p(3)*y(2)^2;
        B2=p(8)*B1;
        F=[A1+B1;A2+B2];        
%         F=[p(1)*y(1) - (p(4)*y(1)*y(2)*(exp(-p(5)*y(1)) - 1))/(p(6) + p(7)*y(2))
%            p(2)*y(2) + p(3)*y(2)^2 - (p(4)*p(8)*y(1)*y(2)*(exp(-p(5)*y(1)) - 1))/ (p(6) + p(7)*y(2))];
% 
     
    case 4
%         dy1dt: '(p1*y1)/(p2 + y1) + (p4*y1*y2)/(p5*y1 + 1)'
%         dy2dt: 'p3*y2 + (p4*p6*y1*y2^2)/((p7 + p8*y2)*(p5*y1 + 1))'

        A1=p(1)*y(1)/(p(2) + y(1));
        B1=p(4)*y(1)*y(2)/(p(5)*y(1) + 1);
        A2=p(3)*y(2);
        B2=p(6)*y(2)/(p(7) + p(8)*y(2))*B1;
        F=[A1+B1;A2+B2];  

%         F=[(p(1)*y(1))/(p(2) + y(1)) + (p(4)*y(1)*y(2))/(p(5)*y(1) + 1)
%            p(3)*y(2) + (p(4)*p(6)*y(1)*y(2)^2)/((p(7) + p(8)*y(2))*(p(5)*y(1) + 1))];

end

end