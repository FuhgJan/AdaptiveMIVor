function [y] = test2_fun(xx)
%% Test 2 in paper
x1 = xx(1);
x2 = xx(2);

fact1 = (sin(x1^2-x2^2))^2 - 0.5;
fact2 = (1 + 0.001*(x1^2+x2^2))^2;

y = (fact1/fact2) - 0.25;

end


