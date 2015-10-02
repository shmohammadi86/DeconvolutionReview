function [estimatedC, trace] = deconvoluteDataset_CVX_final(M, G, varargin)
    n = size(M, 1);
    p = size(M, 2);
    q = size(G, 2);
    
    trace.opt = -1; % For Huber/SVR, it holds the optimal hyper-parameter for the loss function

    args = inputParser;    
    args.addParamValue('loss_fun', 'L2', @(x) ischar(x)); % L1, L2, Hinge, or Huber
    args.addParamValue('NN', true, @(x) islogical(x)); % Should we enforce NonNegativity constraint explicitly in the optimization problem?
    args.addParamValue('STO', false, @(x) islogical(x)); % Should we enforce Sum-To-One constraint explicitly in the optimization problem?
    args.addParamValue('alpha', 0, @(x) isscalar(x) & x >= 0); % relaxation margin for adaptive filtering of Similar Cell Quantity (SCQ) violating features (0 = no margin, otherwise give 1+alpha multiplicative margin)
    args.addParamValue('SCQ_filter_type', 'none', @(x) ischar(x));  % none, G, M, or both
            
    args.addParamValue('regularize', 'none', @(x) ischar(x)); %Type of regularization: none, L2, L1, elastic
    args.addParamValue('lambda1', 0, @(x) isscalar(x) & x >= 0); % weight of L1 Reg     
    args.addParamValue('lambda2', 0, @(x) isscalar(x) & x >= 0); % weight of L2 Reg  

    args.addParamValue('W', speye(n), @(x) ismatrix(x) & nnz(G < 0) == 0); % Weight matrix for to normalize errors (r_i) 
    
    args.addParamValue('knownC', ones(q, p) ./ (p*q), @(x) ismatrix(x)); % True C for all cell-types*samples. It is used for grid-search in Huber/SVR
    
    args.parse(varargin{:});        
    params = args.Results;

    Ref_max = max(G, [], 2);
    Ref_min = min(G, [], 2);
    
    tic;
    switch params.loss_fun
         case 'L2'
            cvx_clear;              
            estimatedC = zeros(q, p);
            for sample_id = 1:p
                m = M(:, sample_id); 
                switch(params.SCQ_filter_type)
                    case 'both'
                        extreme_sample_mask = (1+params.alpha)*Ref_max < m; 
                        extreme_ref_mask = (1+params.alpha)*m < Ref_min; 
                        bad_feature_mask = extreme_sample_mask | extreme_ref_mask;
                    case 'M'
                        bad_feature_mask = (1+params.alpha)*Ref_max < m; 
                    case 'G'
                        bad_feature_mask = (1+params.alpha)*m < Ref_min;
                    case 'none'
                        bad_feature_mask = false(n, 1);
                end
                if(nnz(~bad_feature_mask) ==  0)
                    error('There is no viable feature for sample %s', sample_id);
                else                
                    cvx_begin quiet
                        variable c(q)
                        switch params.regularize
                            case 'none'
                                minimize( norm(params.W*(G*c- m), 2) )
                            case 'L1'
                                minimize( norm(params.W*(G*c- m), 2) + params.lambda1*norm(c, 1)) 
                            case 'L2'
                                minimize( norm(params.W*(G*c- m), 2) + params.lambda2*norm(c, 2))                                                        
                        end
                        subject to
                            if(params.NN)
                                0 <= c;
                            end
                            if(params.STO)                        
                                ones(1, q)*c == 1;
                            end
                    cvx_end
                    estimatedC(:, sample_id) = c;
                end
            end 

         case 'L1'
            cvx_clear;              
            estimatedC = zeros(q, p);
            for sample_id = 1:p
                m = M(:, sample_id); 
                switch(params.SCQ_filter_type)
                    case 'both'
                        extreme_sample_mask = (1+params.alpha)*Ref_max < m; 
                        extreme_ref_mask = (1+params.alpha)*m < Ref_min; 
                        bad_feature_mask = extreme_sample_mask | extreme_ref_mask;
                    case 'M'
                        bad_feature_mask = (1+params.alpha)*Ref_max < m; 
                    case 'G'
                        bad_feature_mask = (1+params.alpha)*m < Ref_min;
                    case 'none'
                        bad_feature_mask = false(n, 1);
                end
                if(nnz(~bad_feature_mask) ==  0)
                    error('There is no viable feature for sample %s', sample_id);
                else                
                    cvx_begin quiet
                        variable c(q)
                        switch params.regularize
                            case 'none'
                                minimize( norm(params.W*(G*c- m), 1) )
                            case 'L1'
                                minimize( norm(params.W*(G*c- m), 1) + params.lambda1*norm(c, 1)) 
                            case 'L2'
                                minimize( norm(params.W*(G*c- m), 1) + params.lambda2*norm(c, 2))                                                        
                        end
                        subject to
                            if(params.NN)
                                0 <= c;
                            end
                            if(params.STO)                        
                                ones(1, q)*c == 1;
                            end
                    cvx_end
                    estimatedC(:, sample_id) = c;
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
                    m = M(:, sample_id); 
                    cvx_begin quiet
                    variable c(q)
                    minimize( sum(huber(G*c- m, half_length)) )
                    subject to
                        if(params.NN)
                            0 <= c;
                        end
                        if(params.STO)                        
                            ones(1, q)*c == 1;
                        end
                    cvx_end                
                    estimatedC(:, sample_id) = c;
                end
                estimatedC(estimatedC < 0) = 0; 
                estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);  

                trace.C{k} = estimatedC;
                trace.RMSD_list_estim(k) = sqrt(mean( ( M(:) - reshape((G*estimatedC), numel(M), 1) ).^2 ));
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
                    m = M(:, sample_id); 
                    cvx_begin quiet
                        variables c(q, 1) zeta_minus(n, 1) zeta_plus(n, 1)
                        minimize (  sum(zeta_plus + zeta_minus) + lambda*0.5*(c'*c)  )
                        subject to
                            (m - G*c)  <= (zeta_plus + loss_epsilon)
                            -(zeta_minus + loss_epsilon) <= (m - G*c)
                            0 <= zeta_minus
                            0 <= zeta_plus   
                            if(params.NN)
                                0 <= c;
                            end
                            if(params.STO)                        
                                ones(1, q)*c == 1;
                            end                            
                    cvx_end                
                    estimatedC(:, sample_id) = c;
                end
                estimatedC(estimatedC < 0) = 0; 
                estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);  

                trace.C{k} = estimatedC;
                trace.RMSD_list_estim(k) = sqrt(mean( ( M(:) - reshape((G*estimatedC), numel(M), 1) ).^2 ));
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

