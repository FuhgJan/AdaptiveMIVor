function [ y, lb, ub ,x, M ] = Test2()


n = 2;
lb = 0*ones(n,1);        % lower bound
ub = 2*ones(n,1); 


% Initial samples
x = scaled_TPLHD(5,lb,ub);  

M = @(xx) test2_fun(xx);

y = zeros(size(x,1),1);
for i=1:size(x,1)
    y(i,1) = M(x(i,:));
end



end

function [y] = test2_fun(xx)
x1 = xx(1);
x2 = xx(2);

fact1 = (sin(x1^2-x2^2))^2 - 0.5;
fact2 = (1 + 0.001*(x1^2+x2^2))^2;

y = (fact1/fact2) - 0.25;

end


