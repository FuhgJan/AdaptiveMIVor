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
global n_of_Variables


%addpath('help_functions')
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

initNumber = size(X,1);

data_errors = Error_Saver();



no_test_points = 5000 * n_of_Variables;
test_points = lhs_scaled(no_test_points,A_sampling(1,:),A_sampling(2,:));

for i=1:no_test_points
    test_points_response(i) = M(test_points(i,:));
    if test_points_response(i) >= class_limit
        classificationVec(i) = 1;
    else
        classificationVec(i) = 0;
    end
end


test_points = scale_vector_to_unity(lb, ub, test_points);
y = test_points_response;

shouldBeClass1 = numel(find(y >= class_limit));
shouldBeClass2 = numel(find(y < class_limit));
close all;

if n_of_Variables ==2
    xd = (max(test_points(:,1))-min(test_points(:,1)))/200;
    yd =(max(test_points(:,2))-min(test_points(:,2)))/200;
    [xq,yq] =        meshgrid(min(test_points(:,1)):xd:max(test_points(:,1)),min(test_points(:,2)):yd:max(test_points(:,2)));
    
    mymap = [0.8 0.8 0.8;
        1 0 0];
    
    %         drawnow
end

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
    
    NumberCorrectClass1 = 0;
    NumberCorrectClass2 = 0;
    for i=1:no_test_points
        
        [metamodel_response(i),~] = metamodel.construct_modell(test_points(i,:));
        if metamodel_response(i) >= class_limit
            classificationVec_pred(i) = 1;
        else
            classificationVec_pred(i) = 0;
        end
        
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
    
    
    
    
    
    if n_of_Variables ==2
        currentSize = size(metamodel.X,1);
        vqreal =griddata(test_points(:,1),test_points(:,2),classificationVec,xq,yq);
        h= figure(2);
        subplot(2,1,1)
        s = surf(xq,yq,vqreal); hold on;
        colormap(mymap);
        title('Target')
        xlabel('$x_{1}$','Interpreter','Latex')
        ylabel('$x_{2}$','Interpreter','Latex')
        set(gca,'FontSize',16)
        s.EdgeColor ='none';
        view(2);
        
        
        vqpred = griddata(test_points(:,1),test_points(:,2),classificationVec_pred,xq,yq);
        subplot(2,1,2)
        
        s = surf(xq,yq,vqpred); hold on;
        scatter3(metamodel.X(1:end-1,1),metamodel.X(1:end-1,2), 100*ones(size(metamodel.X(1:end-1,1),1),1),60,'b','filled'); hold on;
        scatter3(metamodel.X(end,1),metamodel.X(end,2), 100,100,'k','filled'); hold off;
        colormap(mymap);
        title('Prediction and Samples')
        xlabel('$x_{1}$','Interpreter','Latex')
        ylabel('$x_{2}$','Interpreter','Latex')
        set(gca,'FontSize',16)
        set(gcf, 'Position',[850 1000 500 500]);
        s.EdgeColor = 'none';
        view(2);
        drawnow
        
        
    end
    
    stored_metamodels{iter} = metamodel;
    single_errors = {adaptive_method,data_errors};
    
    
    
end


end
