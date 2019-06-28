clear all;
addpath('help_functions')
addpath('response_functions')

% Define exploration rate and decrease factor
global exploration;
global decrease_factor;

exploration = 0.4;
decrease_factor = 1.1;


% Number of initial samples
n = 10;

% Sampling bounds
lb = [0; 0]; %lower bound
ub = [2;2];  %upper bound

A_sampling =[lb';ub'];

% Create initial samples
x = lhs_scaled(n,lb,ub);

% Response function: FEM, FD, ...
M = @(xx) test2_fun( xx );

% Create responses
for i=1:size(x,1)
    y(i,1) = M(x(i,:));
end

% Limit to distinguish classes in 2d
class_limit = 0.0;

% Scale parameter samples to unity
x = scale_vector_to_unity(lb, ub, x);

% Number of adaptive samples
number_of_adaptive_iterations = 3;

% Number of repetitive iterations
max_iteration = 2;


[output] = chooseSamplingMethod(M,x, y, A_sampling, number_of_adaptive_iterations,max_iteration,class_limit);

