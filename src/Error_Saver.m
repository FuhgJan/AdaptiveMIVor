classdef Error_Saver
    properties
        error_data;
    end
    
    methods
        function obj = Error_Saver()
            obj.error_data =  table;
            
        end
        
        function obj = update(obj,it,X,Y,percentAbove, percentBelow)
            MAE_val = MeanAE(X,Y);
            RMAE_val = RMAE(X,Y);
            RMSE_val = RMSE(X,Y);
            R_sq_val = R_sq(X,Y);
            i= size(obj.error_data,1);
            obj.error_data(i+1,:)= {it ,MAE_val , RMAE_val, RMSE_val, R_sq_val, percentAbove, percentBelow };
            if i==0
                obj.error_data.Properties.VariableNames = {'Iterator','MAE','RMAE','RMSE', 'R_sq', 'percentAbove', 'percentBelow'};
            end
        end
        
        
        function plot_data(obj)
            data = obj.error_data.Variables;
            Iterator = data(:,1);
            MeanAE = data(:,2);
            RMAE = data(:,3);
            RMSE = data(:,4);
            R_sq = data(:,5);
            
            figure
            plot(Iterator, MeanAE, 'LineWidth', 2.0); hold on;
            plot(Iterator, RMAE, 'LineWidth', 2.0); hold on;
            plot(Iterator, RMSE, 'LineWidth', 2.0); hold on;
            plot(Iterator, R_sq, 'LineWidth', 2.0); hold off;
            legend('MAE','MeanAE','RMSE', 'R^2');
            xlabel('Iterations');
            
        end
        
    end
    
    
end

function error_val = standard(X)
% root mean squared error
m = numel(X);
mean_response = (1/m)*(sum(X));
error_val = sqrt((1/(m))*(sum((X-mean_response).^2)));
end

function error_val = RMAE(X,Y)
% relative maximum absolute error
error_val = max(abs(X-Y))/standard(X);
end


function error_val = MeanAE(X,Y)
% mean absolute error
m= numel(X);
error_val = (1/m)*(sum(abs(X-Y)));
end


function error_val = RMSE(X,Y)
% root mean squared error
m= numel(X);
error_val = sqrt((1/m)*(sum((X-Y).^2)));
end


function error_val = R_sq(X,Y)
% R_sq score
m= numel(X);
mean_response = (1/m)*(sum(X));
MSE_val = (sum((X-Y).^2));

divisor = sum((X-mean_response).^2);

error_val = 1-(MSE_val/divisor);
end