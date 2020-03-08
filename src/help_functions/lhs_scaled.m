function [X_scaled]=lhs_scaled(n,lb,ub)
% Create n samples with Matlab's Latin hypercube method within bounds lb and ub 
p=length(lb);
[M,N]=size(lb);
if M<N
    lb=lb';
end    
[M,N]=size(ub);
if M<N
    ub=ub';
end
slope=ub-lb;
offset=lb;
SLOPE=ones(n,p);
OFFSET=ones(n,p);
for i=1:p
    SLOPE(:,i)=ones(n,1).*slope(i);
    OFFSET(:,i)=ones(n,1).*offset(i);
end
X_normalized = lhsdesign(n,p);
X_scaled=SLOPE.*X_normalized+OFFSET;

end
