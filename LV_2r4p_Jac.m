function J = LV_2r4p_Jac(t, y, p)
%calculation of ODE jacobian

J = [p(1,1)+p(1,2)*y(2),  p(1,2)*y(1);
    p(2,2)*y(2), p(2,1)+p(2,2)*y(1)];
end