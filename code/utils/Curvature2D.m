function [ crv ] = Curvature2D ( x, y )
% If your curve is given by two row vectors, x and y, you can
% approximate its curvature at each point by the reciprocal of the
% radius of a circumscribing triangle with that point, the preceding
% point, and the succeeding point as vertices. The radius of such a
% triangle is one fourth the product of the three sides divided by its
% area. A zero is appended to the beginning and end of the result to
% make it a vector the same length as x and y.
% 
%   The curvature will be positive for curvature to the left and
% negative for curvature to the right as you advance along the curve.
% You can get the approximate maximum curvature using the 'max' and
% 'abs' functions on the 'crv' vector.

    x1 = x(1:end-2); x2 = x(2:end-1); x3 = x(3:end);
    y1 = y(1:end-2); y2 = y(2:end-1); y3 = y(3:end);
    
    a = sqrt((x3-x2).^2+(y3-y2).^2); % a, b, and c are the three sides
    b = sqrt((x1-x3).^2+(y1-y3).^2);
    c = sqrt((x2-x1).^2+(y2-y1).^2);
    A = 1/2*(x1.*y2+x2.*y3+x3.*y1-x1.*y3-x2.*y1-x3.*y2); % The triangle's area
    crv = [0,4*A./(a.*b.*c),0]; % The reciprocal of its circumscribed radius
end

