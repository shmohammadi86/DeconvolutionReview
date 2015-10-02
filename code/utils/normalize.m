function [P, D]=normalize(A, varargin)
% Normalize columns (rows) of A to have p-norm=1. If column-normalization
% is selected (dim = 1), A = P*D, else if row-normalization
% is selected (dim = 2), A = D*P

    params = inputParser;    
    params.addParamValue('dim', 1      ,@(x) isscalar(x) & x > 0); 
    params.addParamValue('pnorm'         , 1     ,@(x) isscalar(x) & x > 0); 
    params.parse(varargin{:});
    par = params.Results;
    
    if(par.pnorm == Inf)
        A_norm = max(abs(A), par.dim);
    else
        A_norm=sum(abs(A).^par.pnorm, par.dim).^(1/par.pnorm);                
    end
    
    
    if par.dim == 1 %column
        P = A*diag(spfun(@(x) 1./x, A_norm));    
    else
        P = diag(spfun(@(x) 1./x, A_norm))*A;
    end

    
    if(size(A_norm, 1) == 1)
        A_norm = A_norm';
    end
    D = spdiags(A_norm, 0, numel(A_norm), numel(A_norm));        
end