function [ X ] = solveNNLS( A, B, method, X0)
%solveNNLS Uses various method to solve min_{X} ||AX - B||_2^2, s.t. 0 <= X 
    
    if ~exist('method', 'var')
        method = 'nnlsm_blockpivot'; % Default to block pivot method
    end
    if ~exist('X0', 'var')
        X0 = ones(size(A, 2), size(B, 2)); 

    end
    
    switch(method)
        case 'nnlsm_blockpivot'
            X = nnlsm_blockpivot(A, B, 0, X0);
        case 'nnlsm_activeset'
            X = nnlsm_activeset(A, B, 0, X0);
        case 'lsqnonnegvect'         
            X = lsqnonnegvect(A, B);
        case 'largennls'
            X = largennls(A, B, X0);
    end
end

