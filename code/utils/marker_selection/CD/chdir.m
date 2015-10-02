function [ b ] = chdir( ctrl,expm, r )
%       This function caclulates the  characteristic direction for a gene 
%       expression dataset.
%  	ctrl: control gene expressoion data
%  	expm: experiment gene expression data
%  	b: return value, a vector of n-components, representing  characteristic
%          direction of that gene expression dataset. n equals the number of 
%          genes in the expression dataset. b is sorted by its components' 
%          absolute values in descending order.
%  	r: regulaized term. A parameter that smooths the covariance matrix 
%          and reduces potential noise in the dataset.
%
%
% For the input matrix rows should represent genes and columns gene expression 
% profiles.
% r is the regulization term ranging [0,1]. b is the characteristic direction.
% ctrl(control) and expm(experiment) matrices should have the same number
% of genes(rows). 
%
% The function is an implementation of the characteristic direction (chdir) 
% algorithm in the paper. No special toolbox is required to run the function.
%
% Author: Qiaonan Duan
% Ma'ayan Lab, Icahn School of Medicine at Mount Sinai
% Jan.8, 2014
%
% Add gene symbols to results. Apr. 4, 2014
% 

% if no input regulization term, set it to 1
if nargin < 4
    r = 1;
end

err1 = MException('EqualCheck:NotEqual', ...
        'Control expression data must have equal number of genes as experiment expression data!');
if size(ctrl,1)~= size(expm,1)
    throw(err1);
end

err2 =  MException('NaNCheck:HasNaN', ...
        'Control expression data and experiment expression data have to be real numbers. NaN was found');

if any(isnan(ctrl(:))) || any(isnan(expm(:)))
    throw(err2)
end


% There should be variance in expression values of each gene. If  
% gene expression values of a gene are constant, it would dramatically
% affect LDA caculation and results in a wrong answer.
constantThreshold = 1e-5;
ctrlConstantGenes = var(ctrl,0,2) < constantThreshold;
expmConstantGenes = var(expm,0,2) < constantThreshold;

if any(ctrlConstantGenes)
    err3 = MException('ConstantCheck:HasConstantRows',...
        sprintf('%s row(s) in control expression data are constant. Consider Removing the row(s).', num2str(find(ctrlConstantGenes))));
    throw(err3);
elseif any(expmConstantGenes)
    err3 = MException('ConstantCheck:HasConstantRows',...
        sprintf('%s row(s) in experiment expression data are constant. Consider Removing the row(s).', strjoin(find(expmConstantGenes),' ')));
    throw(err3);
end
    
% place control gene expression data and experiment gene expression data into
% one matrix
combinedData = [ctrl expm];

% get the number of samples, namely, the total number of replicates in  control 
% and experiment. 
 [~,samplesCount] = size(combinedData);

 % centralize the matrix by the mean of each row, which is a requirement of PCA.
shiftData = combinedData - repmat( mean(combinedData,2), 1, samplesCount );

% the number of output components desired from PCA. We only want to calculate
% the chdir in a subspace that capture most variance in order to save 
% computation workload. The number is set 20 because considering the number of 
% genes usually present in an expression matrix 20 components would  capture 
% most of the variance.
componentsCount = min([samplesCount-1,20]);

% use the nipals PCA algorithm to calculate R, V, and pcvars. nipals algorithm
% has better performance than the algorithm used by Matlab's builtin PCA function.
% R are scores and V are coefficients or loadings. pcvars are the variances 
% captured by each component 
[R,V,pcvars] = nipals(shiftData',componentsCount,1e5,1e-4);


% We only want components that cpature 95% of the total variance or a little above.
% cutIdx is the index of the compoenent, within which the variance is just equal
% to or a little greater than 95% of the total.
cutIdx = find(cumsum(pcvars) > 0.95);
if isempty(cutIdx)
    cutIdx = componentsCount;
else
    cutIdx = cutIdx(1);
end

% slice R and V to only that number of components.
R = R(:,1:cutIdx);
V = V(:,1:cutIdx);

% the difference between experiment mean and control mean.
meanvec = mean(expm,2) - mean(ctrl,2);

% All the following steps calculate shrunkMats. Refer to the paper for detail.
% ShrunkenMats are the covariance matrix that is placed as denominator 
% in LDA formula. Notice the shrunkMats here is in the subspace of those components
% that capture about 95% of total variance.
Dd = R'*R/samplesCount;
Dd = diag(diag(Dd,0));
sigma = mean(diag(Dd,0));
shrunkMats = r*Dd + sigma*(1-r)*eye(size(R,2));

% The LDA formula.
% V/shrunkMats*V' transforms the covariance matrix from the subspace to full space.
% In Matlab, V/shrunkMats is faster and more accurate than V*inv(shrunkMats).
b = V/shrunkMats*V'*meanvec; 

% normlize b to unit vector
b = b/norm(b);
end

