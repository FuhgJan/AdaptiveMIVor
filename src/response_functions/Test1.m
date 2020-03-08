function [ y, lb, ub ,x, M ] = Test1()


n = 2;
lb = 0*ones(n,1);        % lower bound
ub = 2*ones(n,1); 


% Initial samples
x = scaled_TPLHD(5,lb,ub);  

M = @(xx) test1_fun(xx);

y = zeros(size(x,1),1);
for i=1:size(x,1)
    y(i,1) = M(x(i,:));
end



end


function [y] = test1_fun(xx)
x1 = xx(1);
x2 = xx(2);




if x1 < 0.5
fact2 = -(x1^1+(x2-1)^2);    
else
    fact2 = -3*(x1^2+(x2-1)^2)-0.1;
end
y = fact2 + 0.55;

end


