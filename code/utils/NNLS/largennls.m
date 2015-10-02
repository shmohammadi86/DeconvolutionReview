function res=largennls(E,b,init)
% function res=largennls(E,b)
%
% Vectorized version of the nnls matlab function
% E=g, b=data(M)

tol = 10*eps*norm(E,1)*length(E);
[m,n] = size(E);  %dimensions of g - genes(m) X CTs(n)
itmax = 3*n;

if nargin==3
    tol=tol/100;
end

for i=1:size(b,2)   %i runs on b's columns (# of patients)
   
   f=b(:,i);        %f= the ith column of b
   P = zeros(1,n);  %P = row vector of n zeroes (# of CTs)
   Z = 1:n;         %Z = vector of 1:n (# of CTs)
   if nargin==3
       x = init(:,i);
   else
       x = P';      %x = column vector of n zeroes (# of CTs)
   end  
   ZZ=Z;            %ZZ = vector of 1:n (# of CTs)
   w = E'*(f-E*x);
   iter = 0;
   erriter=1;
   
   while any(Z) & any(w(ZZ) > tol) & erriter;
      [wt,t] = max(w(ZZ));
      t = ZZ(t);
      P(1,t) = t;
      Z(t) = 0;
      PP = find(P);
      ZZ = find(Z);
      nzz = size(ZZ);
      EP(1:m,PP) = E(:,PP);
      EP(:,ZZ) = zeros(m,nzz(2));
      % Ici...
      [U,S,V] = svd(EP,0);
      s = diag(S);
      tol = max(m,n) * max(s) * eps;
      r = sum(s > tol);
      s = diag(ones(r,1)./s(1:r));
      X = V(:,1:r)*s*U(:,1:r)';
      z=X*f;
      %z = pinv(EP)*f;
      z(ZZ) = zeros(nzz(2),nzz(1));
      % inner loop to remove elements from the positive set which no longer belong
      while any((z(PP) <= tol)) & erriter;
         iter = iter + 1;
         if iter > itmax;
            disp(['Iteration count is exceeded.', ...
                  '  Try raising the tolerance.'])
            erriter=0;
         end
         QQ = find((z <= tol) & P');
         alpha = min(x(QQ)./(x(QQ) - z(QQ)));
         x = x + alpha*(z - x);
         ij = find(abs(x) < tol & P' ~= 0);
         Z(ij)=ij';
         P(ij)=zeros(1,length(ij));
         PP = find(P);
         ZZ = find(Z);
         nzz = size(ZZ);
         EP(1:m,PP) = E(:,PP);
         EP(:,ZZ) = zeros(m,nzz(2));
         % Ici...
         [U,S,V] = svd(EP,0);
         s = diag(S);
         tol = max(m,n) * max(s) * eps;
         r = sum(s > tol);
         s = diag(ones(r,1)./s(1:r));
         X = V(:,1:r)*s*U(:,1:r)';
         z=X*f;
         %z = pinv(EP)*f;
         z(ZZ) = zeros(nzz(2),nzz(1));
      end
      x = z;
      w = E'*(f-E*x);
   end
   res(:,i) = x;
end

