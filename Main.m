 clc
 clear

%Init(3);
 %Interpolate('@(x1,x2)sin(x1)+cos(x2)','quad',50,[-1.7,1],[0,3.149]);
 %Interpolate('@(x1,x2)sin(x1)','quad',50,[0,1],[0,1]);
 %Interpolate('@(x1,x2)x1.^4+x2.^4-10','lin',50,[-2,2],[-2,2]);
% Interpolate('@OptBVP_FuncDeriv','quad',5,[-1,1],[-3,2],[0,0]);
%[opx,opf,ope,opout]= fmincon(@OptBVP_FuncDeriv,[6,6,6],[],[],[],[],[-10,-10,-10],[10,10,10])
f=@(x)100.*(x(2)-x(1).^2).^2+(1-x(1)).^2;
df_x1=@(x)-400.*x(1).*x(2)+400.*x(1).^3+2.*x(1)-2;
df_x2=@(x)200.*x(2)-200.*x(1).^2;
%df_x1=@(x)0;
%df_x2=@(x)0;
start_point=[0,0];
lower_bound=[-3,-3];
upper_bound=[3,3];
tic
 Interpolate(f,'quar',10,start_point,lower_bound,upper_bound,df_x1,df_x2);
 toc
 tic
 [x,f,e,out]= fmincon(f,start_point,[],[],[],[],lower_bound,upper_bound)
 toc