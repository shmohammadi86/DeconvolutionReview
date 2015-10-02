    clear
    addpath(genpath('./code/'));
    load('Report_org');
    datasetRange = [2:6];
    Loss_functions = {'L2'};    

    
    current_report = cell(numel(Loss_functions)*4+1, 8);
    current_report(1, :) = {'Loss', 'NN', 'STO', 'mAD', 'mAD pval', 'RMSD', 'RMSD pval', 'dt (s)'};        
    for ds = 1:numel(datasetRange)
        datasetName = datasetNames{datasetRange(ds)};
        fprintf('Loading dataset %s ...\n', datasetName);

        
        datasetPaths = setDatasetPaths(datasetName);
        datasetPaths.annotations.features = ''; % Who cares!
        [ M, G, annotations, known_proportions, H ] = loadDataset(datasetPaths);  

        [n, p] = size(M);
        q = size(G, 2);        



        % Computing weight matrix
        cell_types = unique(annotations.pure.Class);
        prior_C = ones(numel(cell_types), 1) ./ numel(cell_types);

        Gene_variance = cell2mat(arrayfun(@(c) var(H(:, ismember(annotations.pure.Class, cell_types{c})), [], 2), 1:numel(cell_types), 'UniformOutput', false));
        Estimated_mixture_variance = max(Gene_variance, [], 2);%sum(bsxfun(@times, (prior_C.^2)', Gene_variance), 2);

        P = normalize(G, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
        I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
        H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
        
%         x = (1:n);        
%         y = (log2(sort(Estimated_mixture_variance)))';
%         Diff =  smooth(diff(y)/range(y), round(n/100));
%         lower_threshold = 2^y(find(Diff<1e-4, 1));
%         upper_threshold = 2^y(find(Diff<1e-4, 1, 'last'));
%         filter_mask = (Estimated_mixture_variance < lower_threshold) | (upper_threshold < Estimated_mixture_variance);
%         
%         M(filter_mask, :) = [];
%         G(filter_mask, :) = [];
%         H(filter_mask, :) = [];
%         Estimated_mixture_variance(filter_mask, :) = [];
%         Actual_mixture_variance(filter_mask, :) = [];
%         n = nnz(~filter_mask);
        
%         W_std = spdiags(1 ./ std(H, [], 2), 0, n, n);

%         [estimatedC, trace] = deconvoluteDataset_CVX_final(M, G, 'loss_fun', 'L2', 'NN', false, 'STO', false);% , 'regularize', 'L2', 'filter_type', 'none');
%         mAD = sum(sum(abs(estimatedC - known_proportions.C)));
% 
%         Delta{ds} = abs(G*known_proportions.C-M) - abs(G*estimatedC-M);
% %         [CC, PVAL] = corr(Delta)
%         Delta{ds}(:, end+1) = mean(Delta{ds}, 2);
%         Delta{ds}(:, end+1) = std(Delta{ds}, [], 2);
% 
%         Delta{ds}(:, end+1) = zscore(mean(M, 2));
%         Delta{ds}(:, end+1) = zscore(mean(G, 2));
%         
%         Delta{ds}(:, end+1) = zscore(std(M, [], 2));
%         Delta{ds}(:, end+1) = zscore(Estimated_mixture_variance);
% 
%         Delta{ds}(:, end+1) = zscore(H);
%     
%         Facctor_Corr{ds} = partialcorr(Delta{ds}(:, p+1:end));
% %         Facctor_Corr{ds} = corr(Delta{ds}(:, p+1), Delta{ds}(:, p+2:end));
%         [~, perm] = sort(mean(Delta{ds}, 2),'descend');
%         Delta{ds} = Delta{ds}(perm, :);
%         
%         Col_corr{ds} = corr(G);
%         dlmwrite(sprintf('%s.tsv', datasetName), Delta{ds});


        filter_mask = (zscore(H) > -1.95);
        M(filter_mask, :) = [];
        G(filter_mask, :) = [];
        H(filter_mask, :) = [];
        Estimated_mixture_variance(filter_mask, :) = [];
        n = nnz(~filter_mask);

        W_std = spdiags(1 ./ sqrt(Estimated_mixture_variance), 0, n, n);
        
%         lambda_threshold = 10*mean(arrayfun(@(sample) norm(W_sqrt*(G*known_proportions.C(:, sample) - M(:, sample)), 2), 1:size(M, 2)));
%         idx = 1;        
        for Loss = 1:numel(Loss_functions)
            loss_fun = Loss_functions{Loss};
            fprintf('\tUsing Loss: %s\n', loss_fun);
            for NN= 0:1 % NN or not
                if NN == 0, NN_type = 'Imp'; else NN_type = 'Exp'; end
                fprintf('\t\tNN Type = %s\n', NN_type);                    
                for STO = 0:1 % STO or not
                    if STO == 0, STO_type = 'Imp'; else STO_type = 'Exp'; end

                    fprintf('\t\t\tSTO Type = %s ...', STO_type);       
                    tic;
%                     [estimatedC, trace] = deconvoluteDataset_CVX_weighted_regularized(M, G, W_sqrt, 'loss_fun', loss_fun, 'NN', logical(NN), 'STO', logical(STO), 'filter_type', 'both', 'regularize', 'L2', 'lambda2', 10^6);
                    [estimatedC, trace] = deconvoluteDataset_CVX_final(M, G, 'W', W_std, 'loss_fun', loss_fun, 'NN', logical(NN), 'STO', logical(STO), 'SCQ_filter_type', 'none');% , 'regularize', 'L2', 'filter_type', 'none');
                    dt = toc;

                    
                    [~, mAD, RMSD]  = evaluateC(known_proportions.C, estimatedC, 0);
%                     fprintf('\t\t\tmAD = %d\n', mAD);

%                     [overall_score, sample_pvals] = evaluateC_pval(known_proportions.C, estimatedC);
%                     fprintf(' Score = %f\n', overall_score);
                    fprintf(' mAD = %.2f\n', mAD);

                    row_idx = (Loss-1)*4 + NN*2 + STO + 2;                    
                    current_report{row_idx, 1} = loss_fun;
                    current_report{row_idx, 2} = NN_type; 
                    current_report{row_idx, 3} = STO_type; 

                    current_report{row_idx, 4} = mAD;
%                     current_report{row_idx, 5} = decon_pvals.mAD_pval;
                    
                    current_report{row_idx, 6} = RMSD;
%                     current_report{row_idx, 7} = decon_pvals.RMSD_pval;

                    current_report{row_idx, 8} = dt;     
%                     Score(ds, idx) = mAD; idx = idx +1;
                end
            end 
%             idx = 1;
        end  
        Report{ds} = cell2dataset(current_report);
% 
%         dlmcell(sprintf('output/experiment1/%s_%s.tsv', datasetNames{ds}, datestr(now, 'dd-mmm_THHMM')), current_report);
    end
    
    
    Label = cell(4*numel(Loss_functions), 1);
    idx = 1;
    for Loss = 1:numel(Loss_functions)
        loss_fun = Loss_functions{Loss};
        for NN= 0:1 % NN or not
            if NN == 0, NN_type = 'Imp'; else NN_type = 'Exp'; end
            for STO = 0:1 % STO or not
                if STO == 0, STO_type = 'Imp'; else STO_type = 'Exp'; end
                Label{idx} = sprintf('%s (NN=%s, STO=%s)', loss_fun, NN_type, STO_type);
                idx = idx + 1;
            end
        end
    end
    
%     Data = cell2mat(cellfun(@(X) X.mAD(1:4*numel(Loss_functions)), Report_org(datasetRange), 'UniformOutput', false)) - cell2mat(cellfun(@(X) X.mAD(1:4*numel(Loss_functions)), Report, 'UniformOutput', false));
    Data = cell2mat(cellfun(@(X) X.mAD(1:4*numel(Loss_functions)), Report, 'UniformOutput', false));
    Colors = cbrewer('qual', 'Set1', numel(Label));  
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);            
    hArray = bar(Data');
    for i = 1:numel(Label)
        set(hArray(i),'EdgeColor', Colors(i, :), 'EdgeColor', 'k');        
    end
%     set(gca, 'xticklabel', datasetNames(datasetRange),'FontSize', 17, 'FontWeight','bold');
    set(gca,'FontSize',16, 'FontWeight','bold');        
    xticklabel_rotate(1:numel(datasetRange), 90, datasetNames(datasetRange),'FontSize',18, 'FontWeight','bold', 'interpreter', 'none');    
    ylabel('Quality score (Higher the Better)','FontSize', 18, 'FontWeight','bold')        
    legend(Label, 'interpreter', 'none', 'Location', 'EastOutside','FontSize', 18, 'FontWeight','bold');  
%     
%     export_fig('quality_scorE_org.eps', '-eps', '-transparent', '-painters');    
%     close;  

%%
%     % Omega: Multiplicative perturbation matrix
%     [~, worst_m] = max(mean(abs(known_proportions.C - estimatedC)));
%     m_vec = M(:, worst_m); n = numel(m_vec);
%     q = size(G, 2);
%     
%     upper_bound = 1e2; lower_bound = 1e-2;
%     cvx_begin quiet
%         variable Omega(n, q)
%         minimize( norm(G*estimatedC- m_vec, 1) + (1 ) )
%         subject to
% 
%     cvx_end

% cell_types = unique(annotations.pure.Class);
% prior_C = ones(numel(cell_types), 1) ./ numel(cell_types);
% 
% Gene_variance = cell2mat(arrayfun(@(c) var(H(:, ismember(annotations.pure.Class, cell_types{c})), [], 2), 1:numel(cell_types), 'UniformOutput', false));
% Var_of_Var = var(Gene_variance, [], 2);
% ksdensity(log10(sqrt(Var_of_Var)));
% 
% Estimated_mixture_variance = sum(bsxfun(@times, (prior_C.^2)', Gene_variance), 2);
% ksdensity(log10(sqrt(Estimated_mixture_variance)));
% 
% [n, q] = size(G);
% W = spdiags(1 ./ Estimated_mixture_variance, 0, n, n);

% I_q = eye(q);
% lambda = 10^6;
% shifted_product = (G'*W*G + lambda*I_q);
% D = shifted_product\G'*W;
% for i = 1:size(M, 2)
%     m = M(:, i);
%     c = D*m;
%     estimatedC(:, i) = c;
% end
% 
% estimatedC(estimatedC < 0) = 0; 
%  estimatedC = normalize(estimatedC, 'pnorm', 1, 'dim', 1);
