function results = run_benchmark(datasetNames, lossFunctions, outputFile)
% RUN_BENCHMARK  Reproduce the survey's loss x constraint comparison.
%
%   results = run_benchmark()
%   results = run_benchmark(datasetNames)
%   results = run_benchmark(datasetNames, lossFunctions)
%   results = run_benchmark(datasetNames, lossFunctions, outputFile)
%
% For every dataset, every loss function, and each of the four combinations of
% implicit/explicit non-negativity (NN) and sum-to-one (STO) constraints, this
% deconvolutes the mixtures and scores the result against the known proportions.
%
% Inputs (all optional):
%   datasetNames   cell array of names understood by setDatasetPaths.
%                  Default: every dataset with bundled ground truth.
%   lossFunctions  cell array from {'L2', 'L1', 'Huber', 'Hinge'}.
%                  Default: all four.
%   outputFile     path to write the results table as TSV. Default: none.
%
% Output:
%   results   table with one row per (dataset, loss, NN, STO) combination and
%             columns for KL divergence, mean absolute deviation, RMSD, and
%             wall-clock time.
%
% Requires CVX on the path; run cvx_setup once before calling this.
%
% Example:
%   addpath(genpath('code'))
%   r = run_benchmark({'CellLines'}, {'L2', 'L1'}, 'results.tsv');

    if nargin < 1 || isempty(datasetNames)
        datasetNames = {'LiverBrainLung', 'BreastBlood', 'CellLines', ...
                        'RatBrain', 'Retina', 'PERT_Cultured', ...
                        'PERT_Uncultured'};
    end
    if nargin < 2 || isempty(lossFunctions)
        lossFunctions = {'L2', 'L1', 'Huber', 'Hinge'};
    end
    if nargin < 3
        outputFile = '';
    end
    if ischar(datasetNames),   datasetNames   = {datasetNames};   end
    if ischar(lossFunctions),  lossFunctions  = {lossFunctions};  end

    if exist('cvx_begin', 'file') ~= 2
        error('run_benchmark:noCVX', ...
              ['CVX was not found on the path. Install it from ' ...
               'https://cvxr.com/cvx/ and run cvx_setup.']);
    end

    rows = {};

    for d = 1:numel(datasetNames)
        name = datasetNames{d};
        fprintf('=== %s ===\n', name);

        paths = setDatasetPaths(name);
        paths.annotations.features = '';   % not used by this experiment
        [M, G, annotations, known_proportions, H] = loadDataset(paths);

        % Weight each residual by the inverse standard deviation of the gene's
        % expression across cell types, as in the survey's weighted variants.
        cellTypes     = unique(annotations.pure.Class);
        geneVariance  = cell2mat(arrayfun( ...
            @(c) var(H(:, ismember(annotations.pure.Class, cellTypes{c})), [], 2), ...
            1:numel(cellTypes), 'UniformOutput', false));
        mixtureVar    = max(geneVariance, [], 2);

        n = size(M, 1);
        W = spdiags(1 ./ sqrt(max(mixtureVar, eps)), 0, n, n);

        for l = 1:numel(lossFunctions)
            loss = lossFunctions{l};
            for NN = [false true]
                for STO = [false true]
                    fprintf('  %-6s NN=%d STO=%d ... ', loss, NN, STO);
                    t0 = tic;
                    estimatedC = deconvoluteDataset_CVX_final( ...
                        M, G, 'W', W, 'loss_fun', loss, ...
                        'NN', NN, 'STO', STO, 'SCQ_filter_type', 'none');
                    dt = toc(t0);

                    [KL, mAD, RMSD] = evaluateC(known_proportions.C, estimatedC, 0);
                    fprintf('mAD = %.3f  RMSD = %.3f  (%.1fs)\n', mAD, RMSD, dt);

                    rows(end+1, :) = {name, loss, NN, STO, KL, mAD, RMSD, dt}; %#ok<AGROW>
                end
            end
        end
    end

    results = cell2table(rows, 'VariableNames', ...
        {'Dataset', 'Loss', 'NN', 'STO', 'KL', 'mAD', 'RMSD', 'Seconds'});

    if ~isempty(outputFile)
        writetable(results, outputFile, 'FileType', 'text', 'Delimiter', '\t');
        fprintf('\nWrote %s\n', outputFile);
    end
end
