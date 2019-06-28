function x_opti = optimizationTools(fun,strategy,AA,b,Aeq,beq,lb,ub,nonlcon)
% General function that connects to different matlab optimization schemes
n = numel(lb);
if strcmp(strategy,'fmincon')
    options = optimoptions('fmincon','Display','none');
    addpath('help_functions')
    x0 = scale_rand(1,lb,ub);
    x_opti = fmincon(fun,x0,AA,b,Aeq,beq,lb,ub,nonlcon,options);
elseif strcmp(strategy,'PSO')
    options = optimoptions('particleswarm','SwarmSize',800*n,'Display','off');
    options.HybridFcn = @fmincon;
    [x_opti] = particleswarm(fun,n,lb,ub,options);
elseif strcmp(strategy,'GA')
    options = gaoptimset('PopulationSize',10000*n,'Display','off');
    %      options.HybridFcn = @fmincon;
    x_opti = ga(fun,n,AA,b,Aeq,beq,lb,ub,nonlcon,options);
elseif strcmp(strategy,'AN')
    
    options = optimoptions('simulannealbnd','Display','none','FunctionTolerance', 1e-08);
    exitflag = 0;
    while exitflag~= 1
        x0 = (lb+ub)/2;
        [x_opti, ~, exitflag, ~]= simulannealbnd(fun,x0,lb,ub,options);
        ub = ub/1.1;
    end
end
end

