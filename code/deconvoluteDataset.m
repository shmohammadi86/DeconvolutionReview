function [estimatedC] = deconvoluteDataset(M, G, method)
    gene_count = size(M, 1);
    sample_count = size(M, 2);
    CT_count = size(G, 2);
    
    l = 0; u = 1;
    fprintf('Deconvolution using %s method ...\n', method);

    lsqlin_options = optimset('LargeScale','off', 'Display', 'off');        
    switch method
        case 'LS'
            estimatedC = G\M;


         case 'NNLS'
            estimatedC = zeros(CT_count, size(M, 2));
            for sample_id=1:sample_count
                estimatedC(:, sample_id) = lsqnonneg(G, M(:, sample_id));
            end 

         case 'QP'
            estimatedC = zeros(CT_count, size(M, 2));
            for sample_id=1:sample_count
                estimatedC(:, sample_id) = lsqlin(G, M(:, sample_id), [], [], ones(1, CT_count), 1, zeros(CT_count, 1), [], zeros(CT_count, 1), lsqlin_options);
            end 

        case 'nnlsm_blockpivot'
            X0 = ones(CT_count, sample_count);
            X0 = normalize(X0, 'pnorm', 1, 'dim', 1);                    
            X = nnlsm_blockpivot(G, M, 0, X0);
            estimatedC = X;
            
        case 'nnlsm_activeset'
            X0 = ones(CT_count, sample_count);
            X0 = normalize(X0, 'pnorm', 1, 'dim', 1);                    
            X = nnlsm_activeset(G, M, 0, 0, X0);
            estimatedC = X;
            
        case 'largennls'
            X0 = ones(CT_count, sample_count);
            X0 = normalize(X0, 'pnorm', 1, 'dim', 1);                    
            X = largennls(G, M, X0);
            estimatedC = X;
            
        case 'lsqnonnegvect'         
            estimatedC = lsqnonnegvect(G, M);

        case 'SVR'
            C = 1;
            M_norm = M;
            L_norm = G;
            estimatedC = zeros(CT_count, size(M, 2));
            for sample_id=1:sample_count
                m = M_norm(:, sample_id);
                best_rmse = inf; best_w = zeros(CT_count, 1);
                for nu = 0.25:0.25:0.75
                    tic
                    model = svmtrain(m, L_norm, sprintf('-s 4 -t 0 -c %f -n %f -h 0 -m 1024 -e 0.01', C, nu)); % -s: (nu-SVR) -t: linear kernel, -c: sets the parameter C (default 1), -m: 1G memory, -h: no shrinkage heuristic, -e: tolerance=1e-3
                    toc
                    w = (model.sv_coef' * full(model.SVs))'; % based on how the dual form of the SVM optimisation
                    w(w < 0) = 0;
                    w = w ./ sum(w);
                    rmse = sqrt( mean( (L_norm*w - m).^2 ) );
                    if(rmse < best_rmse)
                        best_rmse = rmse;
                        best_w = w;
                    end
                end
                estimatedC(:, sample_id) = best_w;
            end 

        case 'SVR_Z'
            C = 1;
            M_norm = zscore(M);
            L_norm = zscore(G);
            estimatedC = zeros(CT_count, size(M, 2));
            for sample_id=1:sample_count
                m = M_norm(:, sample_id);
                best_rmse = inf; best_w = zeros(CT_count, 1);
                for nu = 0.25:0.25:0.75
                    tic
                    SVR_cmd = sprintf('-s 4 -t 0 -c %f -n %f -h 0 -m 1024 -e 0.01', C, nu); % -s: (nu-SVR) -t: linear kernel, -c: sets the parameter C (default 1), -m: 1G memory, -h: no shrinkage heuristic, -e: tolerance=1e-3
                    model = svmtrain(m, L_norm, SVR_cmd);
                    toc
                    w = (model.sv_coef' * full(model.SVs))'; % based on how the dual form of the SVM optimisation
                    w(w < 0) = 0;
                    w = w ./ sum(w);
                    rmse = sqrt( mean( (L_norm*w - m).^2 ) )
                    if(rmse < best_rmse)
                        best_rmse = rmse;
                        best_w = w;
                    end
                end
                estimatedC(:, sample_id) = best_w;
            end 

        case 'SVR_YALMIP'
            w = sdpvar(CT_count, 1);
            zeta_plus = sdpvar(gene_count, 1);
            zeta_minus = sdpvar(gene_count, 1);

            loss_epsilon = 0; %0.013504;
            lambda = 100;

            M_norm = zscore(M);
            L_norm = zscore(G);
            estimatedC = zeros(CT_count, size(M, 2));
            for sample_id=1:sample_count
                m_norm = M_norm(:, sample_id);

                Objective = sum(zeta_plus + zeta_minus) + lambda*0.5*(w'*w);
                Constraints = [-(zeta_minus + loss_epsilon) <= (m_norm - L_norm*w)  <= (zeta_plus + loss_epsilon), 0 <= zeta_plus, 0 <= zeta_minus, 0 <= w ];

                %options = sdpsettings('verbose', 1, 'solver','cplex','cplex.qpmethod',1);
                options = sdpsettings('verbose', 1);
                sol = optimize(Constraints, Objective, options);
                if sol.problem == 0
                w_val = value(w);
                else
                    display('Hmm, something went wrong!');
                    sol.info
                    yalmiperror(sol.problem)
                end                    
                estimatedC(:, sample_id) = w_val;
            end 

        case 'SVR_CVX'               
            loss_epsilon = 0; %0.013504;
            lambda = 1;

            M_norm = zscore(M);
            L_norm = zscore(G);
            estimatedC = zeros(CT_count, size(M, 2));
            for sample_id=1:sample_count
                m_norm = M_norm(:, sample_id);

                cvx_begin quiet
                    variables w(CT_count, 1) zeta_minus(gene_count, 1) zeta_plus(gene_count, 1)
                    minimize (  sum(zeta_plus + zeta_minus) + lambda*0.5*(w'*w)  )
                    subject to
                        (m_norm - L_norm*w)  <= (zeta_plus + loss_epsilon)
                        -(zeta_minus + loss_epsilon) <= (m_norm - L_norm*w)
                        0 <= zeta_minus
                        0 <= zeta_plus   
                        0 <= w
                cvx_end
                estimatedC(:, sample_id) = w;
            end 

        otherwise
            error('Unknown deconvolution method %s\n', method);
    end

    % Ensure the contraints are satisfied. Even in the methods that are supposed to directly solve constrained optimization, there is always a boundary of error.
    estimatedC(estimatedC < l) = l; 
    estimatedC(estimatedC > u) = u; 
    estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);
end

