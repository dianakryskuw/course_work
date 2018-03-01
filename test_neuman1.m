function test_neuman1()
clc
xinit=linspace(0,1,30);
yinit=[10 0];
solinit=bvpinit(xinit,yinit);
solinit.x
solinit.y

sol=bvp5c(@odefun,@bcfun,solinit);
figure
plot(sol.x, sol.y(1,:),'r')
xlabel('x') ; ylabel('y')
end

function dydx=odefun(x,y)
     dydx=[y(2)
                 y(1)-1];
end

function res=bcfun(y0,ye)
res=[y0(2)+1
         ye(2)+1];
end