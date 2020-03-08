function [ y, lb, ub ,x, M ] = Test3()


n = 2;
lb = [1; 0]; %lower bound
ub = [2;1.95];  %upper bound
% Initial samples
x = scaled_TPLHD(20,lb,ub);  

M = @(xx) test3_fun(xx);

y = zeros(size(x,1),1);
for i=1:size(x,1)
    y(i,1) = M(x(i,:));
end



end

function [y] = test3_fun(xx)
%% Modified DropWave function
x1 = xx(1);
x2 = xx(2);


if x2 > x1
   y = 1- abs(x1 * x2) ;
else
fact1 = 1+ cos(12*sqrt(x1^2 + x2^2));
fact2 = 0.5*(x1^2 + x2^2)+2;

y = -(fact1/fact2) + 0.05;    
end


end
