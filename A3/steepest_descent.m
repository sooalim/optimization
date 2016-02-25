function [x,f_k, a_cc] = steepest_descent(X, f, tol, N)

% f = @(x, y) 100*(y - x.^2).^2 + (1-x).^2
% Inputs
% X: starting value (initial guess x_0)
% f: function
% tol: tolerance
% N: max iterations
% direction p_k = -grad(f(x_k))
% initial step size a_0 = 1
% step size a_k = grad(f(x_k))*p_k/(p_k^T * Q * p_k)
% x_k+1 = x_k + a_k*p_k
format long
if (nargin < 4) 
    tol = 1e-5;
end %if

k = 1; %initialize count
H = eye(2);
x = X;
MaxIter = N;
a_k = 1;
p_k = Grad(x);

while (norm(a_k * p_k, 2) > tol) && (k < MaxIter)
    
    p_k = -1 * Grad(x);
    p_k = p_k/norm(p_k);
    
    a_k = 1;
    beta = 1e-4;
    a_k = linesearch(a_k, beta, p_k, f, x);
    a_cc(k) = norm(a_k*p_k,2);
    x = x + a_k * p_k;
    f_k(k) = f(x(1), x(2));
    k = k+1;
end %while

end
%end steepest_descent

%backtracking algorithm
function alpha = linesearch(a_k, beta, p_k, f, x)
 in = x + a_k * p_k;
 while feval(f, in(1), in(2)) > (feval(f, x(1), x(2)) + 0.1 * a_k * (p_k' * p_k))
    a_k = beta * a_k;
    in = x + a_k * p_k
 end
 alpha = a_k;
end

function out = Grad(x_k)
f_grad = @(x,y) [2*x - 400*x*(- x^2 + y) - 2, - 200*x^2 + 200*y];
out = f_grad(x_k(1), x_k(2));
end


