classdef OK_model
    % Class for ordinary Kriging
    properties
        auto_correlation_function;
        X;
        Y;
        theta_opti_technique;
        
        m;
        F;
        theta;
        R;
    end
    
    
    methods
        %Constructor
        function obj=OK_model(af,x,y,opti)
            obj.auto_correlation_function = af;
            obj.X = x;  % x = [x1 y1 z1; x2 y2 z2]
            obj.Y = y;
            obj.theta_opti_technique = opti;
            
            obj.m = size(x,1);
            obj.F = ones(size(1:obj.m))';
            obj.theta = optimize_theta(obj);
            
            obj.R = compute_R(obj, obj.theta);
        end
        
        
        
        
        
        function R = compute_R(obj, theta)
            
            
            R =  obj.auto_correlation_function(obj.X,obj.X, theta);
            
            
            %Analytical comparison of regularization methods for Gaussian
            %Processes, Mohammadi
            if sum(isnan(R(:)))
                disp('NAN values')
            end
            
            k = cond(R);
            if k > 10^(12)
                min_lambda = min(eigs(R));
                max_lambda = max(eigs(R));
                
                kmax = 10^10;
                tau_sq = (max_lambda - kmax * min_lambda)/(kmax-1);
                
                R = R+ tau_sq*eye(size(R,1));
            end
            
            
        end
        
        
        function beta_hat = compute_beta_hat(obj,R)
            beta_hat = ((obj.F'*(R\obj.F))\obj.F') * (R\obj.Y);
        end
        
        function sigma_sq_hat = compute_sigma_sq_hat(obj,R,beta_hat)
            sigma_sq_hat = (1/obj.m) *(obj.Y- obj.F*beta_hat)' * (R\(obj.Y- obj.F*beta_hat));
        end
        
        
        function r0 = compute_r0(obj,theta,x0)
            r0= obj.auto_correlation_function(obj.X,x0,theta);
        end
        
        function mu_hat = compute_mu_hat(obj,R,beta_hat,x0,theta)
            r0 = compute_r0(obj,theta,x0);
            mu_hat = beta_hat + r0' * (R\(obj.Y - obj.F*beta_hat));
        end
        
        function sigma_Y_sq_hat = compute_sigma_Y_sq_hat(obj,sigma_sq_hat,x0,theta,R)
            r0 = compute_r0(obj,theta,x0);
            u0 = obj.F' * (R\r0) - 1;
            sigma_Y_sq_hat = sigma_sq_hat * (1 - r0' * (R\r0) + u0 * ((obj.F' * (R\obj.F))\u0));
        end
        
        function theta = optimize_theta(obj)
            
            AA = [];
            b = [];
            Aeq = [];
            beq = [];
            
            n = size(obj.X,2);
            for k=1:n
                iter =1;
                clear distance
                for i=1:obj.m
                    for j=1:obj.m
                        if ~(i == j)
                            distance(iter) = abs(obj.X(i,k) - obj.X(j,k));
                            iter = iter +1;
                        end
                    end
                end
                max_distance = max(distance);
                min_distance = min(distance);
                
                lb(k) = 0.05*min_distance;
                if lb(k) == 0.0
                    lb(k) = 10^(-5);
                end
                ub(k) = 5*max_distance;
            end
            
            fun = @obj.computeMLE;
            theta = optimizationTools(fun,obj.theta_opti_technique,AA,b,Aeq,beq,lb,ub,[]);
            
        end
        
        
        function Psi = computeMLE(obj,theta)
            R_matrix = compute_R(obj, theta);
            beta_hat = compute_beta_hat(obj,R_matrix);
            
            sigma_sq_hat = compute_sigma_sq_hat(obj,R_matrix,beta_hat);
            
            Psi = 0.5 * (obj.m*log(sigma_sq_hat) + log(det(R_matrix)));
            
        end
        
        
        
        
        function [mu_hat,sigma_Y_sq_hat] = construct_modell(obj,x0)
            
            beta_hat = compute_beta_hat(obj,obj.R);
            sigma_sq_hat = compute_sigma_sq_hat(obj,obj.R,beta_hat);
            for i=1:numel(x0(:,1))
                mu_hat(i) = compute_mu_hat(obj,obj.R,beta_hat,x0(i,:),obj.theta);
                sigma_Y_sq_hat(i) = compute_sigma_Y_sq_hat(obj,sigma_sq_hat,x0(i,:),obj.theta,obj.R);
            end
        end
        
        
        %% adaptive sampling methods
        
        function x_new = adaptive_sampling(obj,method,A,class_Limit)
            if strcmp(method,'MiVor')
                x_new = MiVor_function(obj,A,class_Limit);
            elseif strcmp(method,'EIGF')
                x_new = doEIGF(obj,A);
            end
            
            function [S] = SFVCT_S(obj)
                for i=1:obj.m
                    clear distance_min
                    iter = 1;
                    for j=1:obj.m
                        
                        if ~(i==j)
                            distance_min(iter) = norm(obj.X(i,:) - obj.X(j,:));
                            iter = iter +1;
                        end
                    end
                    distance(i) = min(distance_min);
                end
                
                max_distance = max(distance);
                
                S= 0.5 *max_distance;
            end
            
            function [newPoint,points] = findNewPointInHighestVolume(obj,C, indexOfInterest)
                X_without_indexOfInterest = obj.X;
                Y_without_indexOfInterest = obj.Y;
                X_without_indexOfInterest(indexOfInterest,:) = [];
                Y_without_indexOfInterest(indexOfInterest,:) = [];
                [indexNeighbor,~] = knnsearch(X_without_indexOfInterest,obj.X(indexOfInterest,:),'K',2*size(X_without_indexOfInterest,2));
                
                points = C{indexOfInterest,2};
                % Check if one of the closest points is lower Limit
                closestPoints_class2 = 0;
                for kk=1:numel(indexNeighbor)
                    if Y_without_indexOfInterest(indexNeighbor(kk)) < 0
                        class2_pointOfINterest = indexNeighbor(kk);
                        closestPoints_class2 = 1;
                        break;
                    end
                end
                
                if closestPoints_class2
                    disp('Closest point lower limit found!')
                    Distance = inf;
                    for ll=1:size(points,1)
                        d = norm(X_without_indexOfInterest(class2_pointOfINterest) - points(ll,:));
                        if d< Distance
                            Distance = d;
                            newPointTemp = points(ll,:);
                        end
                    end
                else
                    disp('Variance')
                    beta_hat = compute_beta_hat(obj,obj.R);
                    sigma_sq_hat = compute_sigma_sq_hat(obj,obj.R,beta_hat);
                    %                 Distance = inf;
                    sigma_Y_sq_hat_Max = -inf;
                    for ll=1:size(points,1)
                        sigma_Y_sq_hat = compute_sigma_Y_sq_hat(obj,sigma_sq_hat,points(ll,:),obj.theta,obj.R);
                        if sigma_Y_sq_hat > sigma_Y_sq_hat_Max
                            sigma_Y_sq_hat_Max = sigma_Y_sq_hat;
                            newPointTemp = points(ll,:);
                        end
                    end
                end
                newPoint = newPointTemp;
            end
            
            
            
            function x_new = DoMIPT(obj,A)
                addpath('help_functions')
                n = numel(A(1,:));
                lb = A(1,:);
                ub = A(2,:);
                p = lhs_scaled(500 *n* obj.m,lb,ub);
                
                alpha = 0.5;
                dmin = (2*alpha)/ obj.m;
                MIPT_val = -inf;
                
                
                for i=1:size(p,1)
                    val = intersite_proj_th(dmin, obj.X,p(i,:));
                    if val > MIPT_val
                        MIPT_val = val;
                        x_new = p(i,:);
                    end
                end
            end
            
            function x_new = MiVor_function(obj,A, class_Limit)
                addpath('help_functions')
                global exploration;
                global decrease_factor;
                IndexClass1 = find(obj.Y >=class_Limit);
                IndexClass2 = find(obj.Y <class_Limit);
                
                figure(1)
                voronoi(obj.X(:,1),obj.X(:,2)); hold on;
                scatter(obj.X(IndexClass1,1),obj.X(IndexClass1,2),'b'); hold on;
                scatter(obj.X(IndexClass2,1),obj.X(IndexClass2,2),'y'); hold off;
                drawnow
                
                
                r = rand();
                if isempty(IndexClass1) || (r<exploration)
                    disp('Exploration')
                    ST = ['Exploration rate:   ', num2str(exploration)];
                    disp(ST);
                    newPoint = DoMIPT(obj,A);
                    exploration = exploration/decrease_factor;
                else
                    disp('Exploitation')
                    ST = ['Exploration rate:   ', num2str(exploration)];
                    disp(ST);
                    
                    lb = A(1,:);
                    ub = A(2,:);
                    
                    voronoi_matrix = randomVoronoi(obj.X,lb,ub);
                    for I_sample_P=1:obj.m
                        Vol(I_sample_P) = 0.0;
                        for i=2:size(voronoi_matrix{I_sample_P,2},1)
                            %                     Vol(I_sample_P) = Vol(I_sample_P) + (1/(norm(C{I_sample_P,2}(i,:))));
                            Vol(I_sample_P) = Vol(I_sample_P) + 1;
                        end
                    end
                    abs_Vol= sum(Vol);
                    Vol = (1/abs_Vol)*Vol;
                    
                    
                    [~,IndexmaxVol_SubSet] = max(Vol(IndexClass1));
                    indexmaxVol = IndexClass1(IndexmaxVol_SubSet);
                    
                    [newPoint,~] = findNewPointInHighestVolume(obj,voronoi_matrix,indexmaxVol);
                    [~,D] = knnsearch(obj.X,newPoint,'K',obj.m);
                    
                    S = 0.25*SFVCT_S(obj);
                    
                    if numel(find(D< S))>1
                        newPoint = DoMIPT(obj,A);
                    end
                    
                    
                end
                x_new = newPoint;
                
                figure(1)
                voronoi(obj.X(:,1),obj.X(:,2)); hold on;
                scatter(x_new(:,1),x_new(:,2),'r'); hold on;
                scatter(obj.X(IndexClass1,1),obj.X(IndexClass1,2),'b'); hold on;
                scatter(obj.X(IndexClass2,1),obj.X(IndexClass2,2),'y'); hold off;
                drawnow
                
            end
            
            
            
            
            %% EIGF
            function x_new = doEIGF(obj,A)
                AA = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = A(1,:);
                ub = A(2,:);
                strategy = 'AN';
                
                fun = @(x) adaptive_EIGF(obj,x);
                
                x_new = optimizationTools(fun,strategy,AA,b,Aeq,beq,lb,ub,[]); 
                
            end
            
            
            function EIGF_min = adaptive_EIGF(obj,x)

                beta_hat = compute_beta_hat(obj,obj.R);
                sigma_sq_hat = compute_sigma_sq_hat(obj,obj.R,beta_hat);
                
                mu_hat = compute_mu_hat(obj,obj.R,beta_hat,x,obj.theta);
                sigma_Y_sq_hat = compute_sigma_Y_sq_hat(obj,sigma_sq_hat,x,obj.theta,obj.R);
                
                
                
                k = dsearchn(obj.X,x);
                EIGF = (mu_hat - obj.Y(k))^(2) + sigma_Y_sq_hat;
                
                EIGF_min=-EIGF;
            end
            
            
            
        end
    end
end


