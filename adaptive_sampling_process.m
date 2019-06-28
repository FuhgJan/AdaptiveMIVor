function [stored_metamodels,single_errors] = adaptive_sampling_process(metamodel_ini,M,adaptive_method,A_sampling,number_of_adaptive_iterations, class_limit)
% inputs: 
% metamodel_ini - initial metamodel
% M - respone function
% adaptive_method - adaptive sampling technique
% A_sampling - parameter space
% number_of_adaptive_iterations - maximum number of iterations
% class_limit - Sets limit value to distinguish classes in 2d 
%
% outputs: 
% stored_metamodels - stored metamodels
% single_errors - errors of each adaptive step



addpath('help_functions')
iter = 1;
stored_metamodels{1} = metamodel_ini;
metamodel = metamodel_ini;


lb = A_sampling(1,:);
ub = A_sampling(2,:);

lb_unity = zeros(size(lb));
ub_unity = ones(size(ub));
A_sampling_unity =[lb_unity;ub_unity];

Y = metamodel_ini.Y;
X = metamodel_ini.X;


data_errors = Error_Saver();

while (iter < number_of_adaptive_iterations)
    clear y_new x_new
    x_new = metamodel.adaptive_sampling(adaptive_method,A_sampling_unity, class_limit);
    
    ST = ['New found point:   ', num2str(x_new)];
    disp(ST);
    
    x_new_scaled = scale_vector_from_unity(lb,ub,x_new);
    
    
    for ss=1:size(x_new_scaled,1)
        y_new(ss,1) = M(x_new_scaled(ss,:));
    end
    
    Y = [Y; y_new];
    X = [X; x_new];
    
    metamodel = OK_model(metamodel.auto_correlation_function,X,Y,metamodel.theta_opti_technique);
    
    
    iter = iter+1;
    
    %% Errors
    n_of_Variables = size(metamodel.X,2);
    no_test_points = 5000 * n_of_Variables;
    test_points = lhs_scaled(no_test_points,A_sampling(1,:),A_sampling(2,:));
    
    for i=1:no_test_points
        test_points_response(i) = M(test_points(i,:));
    end
    
    
    test_points = scale_vector_to_unity(lb, ub, test_points);
    y = test_points_response;
    
    shouldBeClass1 = numel(find(y >= class_limit));
    shouldBeClass2 = numel(find(y < class_limit));
    NumberCorrectClass1 = 0;
    NumberCorrectClass2 = 0;
    for i=1:no_test_points
        
        [metamodel_response(i),~] = metamodel.construct_modell(test_points(i,:));
        if (test_points_response(i) >= class_limit) && (metamodel_response(i) >= class_limit)
            NumberCorrectClass1 = NumberCorrectClass1+1;
        elseif (test_points_response(i) < class_limit) && (metamodel_response(i) < class_limit)
            NumberCorrectClass2 = NumberCorrectClass2 +1;
        end
    end
    
    PercentAbove = NumberCorrectClass1/shouldBeClass1;
    PercentBelow = NumberCorrectClass2/shouldBeClass2;
    
    data_errors=data_errors.update(metamodel.m,test_points_response , metamodel_response, PercentAbove, PercentBelow );
    ST =[adaptive_method, '   m = ', num2str(size(X,1)), '  Average RMSE value  :     ', num2str(data_errors.error_data.RMSE(end))];
    disp(ST);
    ST =[adaptive_method, '   m = ', num2str(size(X,1)), '  Correct Class 1 in %:     ', num2str(PercentAbove*100)];
    disp(ST);
    ST =[adaptive_method, '   m = ', num2str(size(X,1)), '  Correct Class 2 in %:     ',  num2str(PercentBelow*100)];
    disp(ST);
    
    
    stored_metamodels{iter} = metamodel;
    single_errors = {adaptive_method,data_errors};
    
    
    
end


end
