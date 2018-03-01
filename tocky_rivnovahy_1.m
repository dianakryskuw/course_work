syms p1 p2 p3 p4 p5  p6 y1 y2

s1='p1*y1+p2*y1*y1+p3*y1*y2'
s2='p4*y2+p5*y2*y2+p6*y1*y2'

% s1='p1*y1+p2*y1*y2'
% s2='p3*y2+p4*y2*y2+p5*y1*y2'
sol=solve(s1,s2,'y1','y2')

% p=[p1 p2; p3 p4];
% s1='p(1,1)*y1+p(1,2)*y1*y2'
% s2='p(2,1)*y2+p(2,2)*y1*y2'
% sol=solve(s1,s2,'y1','y2')

disp('sol.y1')
disp(sol.y1)
disp('sol.y2')
disp(sol.y2)

% sol.y1
%       0
%  -p3/p4
% sol.y2
%       0
%  -p1/p2

df=[diff(s1,y1,1) diff(s1,y2,1);
          diff(s2,y1,1) diff(s2,y2,1)]

% df =
% [ p1 + p2*y2,      p2*y1]
% [      p4*y2, p3 + p4*y1]

%  df_2=[p1 + p2*sol.y2(2),      p2*sol.y1(2);
%        p4*sol.y2(2),           p3 + p4*sol.y1(2)]

df_1=subs(df,{y1,y2},{sol.y1(1),sol.y2(1)})
[v_1,eg_1]=eig(df_1) 
pretty(eg_1)

df_2=subs(df,{y1,y2},{sol.y1(2),sol.y2(2)})

% df_2=[p(1,1) + p(1,2)*sol.y2(2),      p(1,2)*sol.y1(2);
%         p(2,2)*sol.y2(2),           p(2,1) + p(2,2)*sol.y1(2)]
[v_2,eg_2]=eig(df_2)    
pretty(eg_2)
%eg_2_n=subs(eg_2,'p1',3)
%eg_2_n=subs(eg_2,'p1','a')
eg_2_n=subs(eg_2,{p1,p2,p3,p4},{3,2,27,4})
%eg_2_n=subs(eg_2,{p1,p2,p3,p4},{'a',2,'b',4})
 
df_3=subs(df,{y1,y2},{sol.y1(3),sol.y2(3)})
[v_3,eg_3]=eig(df_3) 
pretty(eg_3)

df_4=subs(df,{y1,y2},{sol.y1(4),sol.y2(4)})
[v_4,eg_4]=eig(df_4) 
pretty(eg_4)

syms k m r
eg_4_n=subs(eg_4,{p1,p2,p3,p4,p5,p6},{1,-1,-2,-1,-1,3})

eg_4_n=subs(eg_4,{p1,p2,p3,p4,p5,p6},{2.1,-1,-1.5,-1,-1,0.7})

eg_4_n=subs(eg_4,{p1,p2,p3,p4,p5,p6},{2.1,-1,-1.5,-1,-1,r})

eg_4_n=subs(eg_4,{p1,p2,p3,p4,p5,p6},{1,0,-1.5,-1,0,r})
pretty(eg_4_n)
simplify(eg_4_n)
simplify(eg_4_n)