clear all;
addpath(genpath('src'))

% Define exploration rate and decrease factor
global exploration;
global decrease_factor;
global n_of_Variables

exploration = 0.15;
decrease_factor = 1.1;


n_of_Variables = 2;

%% Test function
% Test Case 1
[ y, lb, ub ,x, M ] = Test1();
% Test Case 2
% [ y, lb, ub ,x, M ] = Test2();
% Test Case 3
% [ y, lb, ub ,x, M ] = Test2();




%% Create responses
for i=1:size(x,1)
    y(i,1) = M(x(i,:));
end

% Limit to distinguish classes in 2d
class_limit = 0.0;

% Scale parameter samples to unity
x = scale_vector_to_unity(lb, ub, x);

% Number of adaptive samples
number_of_adaptive_iterations = 150;

% Number of repetitive iterations
max_iteration = 1;


A_sampling =[lb';ub'];
[output] = chooseSamplingMethod(M,x, y, A_sampling, number_of_adaptive_iterations,max_iteration,class_limit);
save('SavedOutput.mat','output')
