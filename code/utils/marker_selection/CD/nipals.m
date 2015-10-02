function [T,P,pcvar] = nipals(X,a,it,tol)

% Nipals algorithm for PCA.
% This function is written largely based on nipals function from R chemometrics package.

	if nargin == 2
		it = 10;
		tol = 1e-4;
	elseif nargin == 3
		tol = 1e-4;
	end

	[obsCount,varCount] = size(X);
	Xh = X - repmat( mean(X,1), obsCount, 1 );
	T = zeros(obsCount,a);
	P = zeros(varCount,a);
	pcvar = zeros(varCount,1);
	varTotal = sum(var(Xh));
	currVar = varTotal;
	nr = 0;

	for h = 1:a
		th = Xh(:,1);
		ende = false;

		while(~ende)
			nr = nr+1;
			ph = Xh'*th/(th'*th);
			ph = ph/norm(ph);
			thnew = Xh*ph/(ph'*ph);
			prec = (thnew-th)'*(thnew-th);
			th = thnew;

			if prec <= tol^2
				ende = true;
			elseif it <= nr
				ende = true;
				disp('Iteration stops without convergence');
			end
		end

		Xh = Xh-th*ph';
		T(:,h) = th;
		P(:,h) = ph;
		oldVar = currVar;
		currVar = sum(var(Xh));
		pcvar(h) = ( oldVar - currVar )/varTotal;
		nr = 0;
	end

				