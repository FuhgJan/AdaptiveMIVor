function [stored_metamodels, final_errors, single_errors] = initial_metamodel(M,A_sampling,adaptive_method, x, y,max_iteration,class_limit)
% inputs: 
% M - response function
% A_sampling - parameter space
% adaptive_method - chosen adaptive method (MiVor here)
% x - samples in parametric space
% y - observations
% max_iteration - number of repetitions
% class_limit - Sets limit value to distinguish classes in 2d 
%
% outputs: 
% stored_metamodels - stored metamodels over adaptive process
% final_errors - final errors after adaptive process
% single_errors - single errors after each adaptive step
n_of_Variables = size(x,2);

%addpath('Autocorrelation_functions')
af = @Matern32_matrix;
metamodel_ini = OK_model(af,x,y,'AN');

if strcmp(adaptive_method,'Initial_error')
    stored_metamodels{1} = metamodel_ini;
    single_errors = 0;
else
    [stored_metamodels,single_errors] = adaptive_sampling_process(metamodel_ini,M,adaptive_method,A_sampling,max_iteration, class_limit);
end
data_errors = Error_Saver();

lb = A_sampling(1,:);
ub = A_sampling(2,:);

no_test_points = 5000 * n_of_Variables;
test_points = lhs_scaled(no_test_points,lb,ub);

for i=1:no_test_points
    test_points_response(i) = M(test_points(i,:));
end

test_points = scale_vector_to_unity(lb, ub, test_points);

%% Errors
shouldBeClass1 = numel(find(test_points_response >= class_limit));
shouldBeClass2 = numel(find(test_points_response < class_limit));
NumberCorrectClass1 = 0;
NumberCorrectClass2 = 0;
for i=1:no_test_points
    [metamodel_response(i),~] = stored_metamodels{numel(stored_metamodels)}.construct_modell(test_points(i,:));
    if (test_points_response(i) >= class_limit) && (metamodel_response(i) >= class_limit)
        NumberCorrectClass1 = NumberCorrectClass1+1;
    elseif (test_points_response(i) < class_limit) && (metamodel_response(i) < class_limit)
        NumberCorrectClass2 = NumberCorrectClass2 +1;
    end
end

PercentAbove = NumberCorrectClass1/shouldBeClass1;
PercentBelow = NumberCorrectClass2/shouldBeClass2;


final_errors=data_errors.update(stored_metamodels{numel(stored_metamodels)}.m,test_points_response , metamodel_response,PercentAbove, PercentBelow);

end


