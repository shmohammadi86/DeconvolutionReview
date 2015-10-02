function [estimatedC, trace] = deconvoluteDataset_CVX_adaptive(Y, X, varargin)
    n = size(Y, 1);
    p = size(Y, 2);
    q = size(X, 2);
    
    trace.opt = -1; % For Huber/SVR, it holds the optimal hyper-parameter for the loss function

    args = inputParser;    
    args.addParamValue('loss_fun', 'L2', @(x) ischar(x)); 
    args.addParamValue('NN', true, @(x) islogical(x));
    args.addParamValue('STO', false, @(x) islogical(x));
    args.addParamValue('knownC', ones(q, p) ./ (p*q), @(x) ismatrix(x)); % True C for all cell-types*samples. It is used for grid-search in Huber/SVR
    args.addParamValue('lambda1', 0, @(x) isscalar(x) & x >= 0); % weight of L1 Reg     
    args.addParamValue('lambda2', 0, @(x) isscalar(x) & x >= 0); % weight of L2 Reg  
    args.addParamValue('alpha', 0, @(x) isscalar(x) & x >= 0); % relaxation for adaptive range filtering
    args.addParamValue('filter_type', 'both', @(x) ischar(x)); 
    
    args.parse(varargin{:});        
    params = args.Results;
    
    Ref_max = max(X, [], 2);
    Ref_min = min(X, [], 2);
    
    tic;
    switch params.loss_fun
         case 'L2'
            cvx_clear;              
            estimatedC = zeros(q, p);
            for sample_id = 1:p
                y = Y(:, sample_id);
                switch(params.filter_type)
                    case 'both'
                        extreme_sample_mask = (1+params.alpha)*Ref_max < y; 
                        extreme_ref_mask = (1+params.alpha)*y < Ref_min; 
                        bad_feature_mask = extreme_sample_mask | extreme_ref_mask;
                    case 'M'
                        bad_feature_mask = (1+params.alpha)*Ref_max < y; 
                    case 'G'
                        bad_feature_mask = (1+params.alpha)*y < Ref_min;  
                 end
                if(nnz(~bad_feature_mask) ==  0)
                    error('There is no viable feature for sample %s', sample_id);
                else
                    cvx_begin quiet
                        variable w(q)
                        minimize( norm(X(~bad_feature_mask, :)*w- y(~bad_feature_mask), 2) )
                        subject to
                            if(params.NN)
                                0 <= w;
                            end
                            if(params.STO)                        
                                ones(1, q)*w == 1;
                            end
                    cvx_end
                    estimatedC(:, sample_id) = w;
                end
            end 

         case 'L1'
            cvx_clear;              
            estimatedC = zeros(q, p);
            for sample_id = 1:p
                y = Y(:, sample_id); 
                extreme_sample_mask = (1+params.alpha)*Ref_max < y;
                extreme_ref_mask = (1+params.alpha)*y < Ref_min;
                filtered.M(sample_id) = nnz(extreme_sample_mask);
                filtered.G = nnz(extreme_ref_mask);                
                bad_feature_mask = extreme_sample_mask | extreme_ref_mask;
                if(nnz(~bad_feature_mask) ==  0)
                    error('There is no viable feature for sample %s', sample_id);
                else
                    cvx_begin quiet
                        variable w(q)
                        minimize( norm(X(~bad_feature_mask, :)*w- y(~bad_feature_mask), 1) )
                        subject to
                            if(params.NN)
                                0 <= w;
                            end
                            if(params.STO)                        
                                ones(1, q)*w == 1;
                            end
                    cvx_end
                    estimatedC(:, sample_id) = w;
                end
            end
            
        case 'Huber'
            best_RMSE = Inf;
            cvx_clear;              
            HalfLengths = logspace(-7, 7, 15);
            trace.RMSD_list_estim = zeros(numel(HalfLengths), 1);
            trace.RMSD_list_actual = zeros(numel(HalfLengths), 1);
            trace.mAD_list = zeros(numel(HalfLengths), 1);
            trace.KL_list = zeros(numel(HalfLengths), 1);
            trace.C = cell(numel(HalfLengths), 1);
            
            estimatedC = zeros(q, p);            
            for k = 1:numel(HalfLengths)
                half_length = HalfLengths(k);
                for sample_id = 1:p
                    y = Y(:, sample_id); 
                    cvx_begin quiet
                    variable w(q)
                    minimize( sum(huber(X*w- y, half_length)) )
                    subject to
                        if(params.NN)
                            0 <= w;
                        end
                        if(params.STO)                        
                            ones(1, q)*w == 1;
                        end
                    cvx_end                
                    estimatedC(:, sample_id) = w;
                end
                estimatedC(estimatedC < 0) = 0; 
                estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);  

                trace.C{k} = estimatedC;
                trace.RMSD_list_estim(k) = sqrt(mean( ( Y(:) - reshape((X*estimatedC), numel(Y), 1) ).^2 ));
                [trace.KL_list(k), trace.mAD_list(k), trace.RMSD_list(k)] = evaluateC(params.knownC, estimatedC, 0);
                
                if(trace.RMSD_list_actual(k) < best_RMSE)
                    trace.opt = half_length;
                    best_RMSE = trace.RMSD_list_actual(k);
                    best_C = estimatedC;
                end
            end
            estimatedC = best_C;
            
        case 'Hinge'   
            % Store trace, and optimal value
            best_RMSE = Inf;
            cvx_clear;              
            lambda = params.lambda2;
            Epsilons = logspace(-7, 7, 15);
            trace.RMSD_list_estim = zeros(numel(Epsilons), 1);
            trace.RMSD_list_actual = zeros(numel(Epsilons), 1);
            trace.mAD_list = zeros(numel(Epsilons), 1);
            trace.KL_list = zeros(numel(Epsilons), 1);
            trace.C = cell(numel(Epsilons), 1);
            
            estimatedC = zeros(q, p);            
            for k = 1:numel(Epsilons)
                loss_epsilon = Epsilons(k);
                for sample_id = 1:p
                    y = Y(:, sample_id); 
                    cvx_begin quiet
                        variables w(q, 1) zeta_minus(n, 1) zeta_plus(n, 1)
                        minimize (  sum(zeta_plus + zeta_minus) + lambda*0.5*(w'*w)  )
                        subject to
                            (y - X*w)  <= (zeta_plus + loss_epsilon)
                            -(zeta_minus + loss_epsilon) <= (y - X*w)
                            0 <= zeta_minus
                            0 <= zeta_plus   
                            if(params.NN)
                                0 <= w;
                            end
                            if(params.STO)                        
                                ones(1, q)*w == 1;
                            end                            
                    cvx_end                
                    estimatedC(:, sample_id) = w;
                end
                estimatedC(estimatedC < 0) = 0; 
                estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);  

                trace.C{k} = estimatedC;
                trace.RMSD_list_estim(k) = sqrt(mean( ( Y(:) - reshape((X*estimatedC), numel(Y), 1) ).^2 ));
                [trace.KL_list(k), trace.mAD_list(k), trace.RMSD_list(k)] = evaluateC(params.knownC, estimatedC, 0);
                
                if(trace.RMSD_list_actual(k) < best_RMSE)
                    trace.opt = loss_epsilon;
                    best_RMSE = trace.RMSD_list_actual(k);
                    best_C = estimatedC;
                end
            end
            estimatedC = best_C;

        otherwise
            error('Unknown loss function %s\n', loss_fun);                        
    end
    estimatedC(estimatedC < 0) = 0; 
    estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);
    
    trace.dt = toc;       
end

