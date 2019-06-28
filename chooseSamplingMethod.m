function [output] = chooseSamplingMethod(M, x, y,A_sampling, number_of_adaptive_iterations,max_iteration, class_limit)
% inputs: 
% M - response function
% x - samples in parametric space
% y - observations
% A_sampling - parameter space
% number_of_adaptive_iterations - maximum number of iterations
% max_iteration - number of repetitions
% class_limit - Sets limit value to distinguish classes in 2d 
%
% outputs: 
% output - adaptive method, stored metamodels, final errors, single error

% MiVor method
adaptive_methods{1} ='MiVor';

% Expected improvement for global fit
adaptive_methods{2} ='EIGF';

% Fit Kriging to initial dataset and obtain error values
adaptive_methods{3} ='Initial_error';

% Choose one or multiple sampling techniques by setting j value. Can be
% parallelized via parfor.
for j=1:1
    close all
    for i=1:max_iteration
        St = ['Computation of ', adaptive_methods{j},'   Iteration number: ', num2str(i), ' of ', num2str(max_iteration)];
        disp(St)
        
        [stored_metamodels{i}, final_errors{i}, single_errors{1,i}] = initial_metamodel(M,A_sampling,adaptive_methods{j}, x, y,number_of_adaptive_iterations, class_limit);
        
    end
    
    output{j,1}= {adaptive_methods{j}, stored_metamodels, final_errors,single_errors};
    
end


end

