function [x,f_k, x_x, x_y, steps, norm_grad] = BFGS(x_0, f, tol, MaxIter, beta, plot)
%[a,f_k,x,y] = BFGS([1.2, 1.2], f, 1e-6,7000, 0.1); 
%(rosenbrock) (196,872 - non normalized)
%[a,f_k,x,y] = BFGS([500,370], f, 1e-6, 7000, 0.1);
%(egg holder, pk = pk/norm(pk)
% f = @(x, y) 100*(y - x.^2).^2 + (1-x).^2 (rosenbrock)
% f = @(x, y) -1*(y+47)*sin(sqrt(abs(x/2+(y+47))))-x*sin(sqrt(abs(x-(y+47)))); % eggholder function
% Inputs
% X: starting value (initial guess x_0)
% f: function
% tol: tolerance
% N: max iterations
% direction p_k = -grad(f(x_k))
% initial step size a_0 = 1
% step size a_k deteremined by backtracking
% x_k+1 = x_k + a_k*p_k
format long

if (nargin < 4) 
    tol = 1e-6;
end %if

k = 1; %initialize count
x=x_0;
% initialize Hessian
% B = eye(2,2);
 B = Hess(x);
 H = inv(B);

%f'(x)
grad1 = Grad(x);

while  (k <= MaxIter) && (norm(grad1,2) > tol) 
    p_k = -1 * H * grad1';
    p_k = p_k/norm(p_k);    
    
    a_k = 1;
    %a_k = linesearch(a_k, beta, p_k, f, x);
    a_k = inexact(a_k, p_k, f, x, grad1);    
    
    s_k = a_k * p_k; % x_(k+1) - x_k
    s_k = s_k'; 
    old_x = x;
    x = x + s_k;
    steps(k) = norm(x - old_x, 2);
    norm_grad(k) = norm(grad1,2);
   
    %{
    if(x(1) < -512 || x(2) < -512 || x(1) > 512 || x(2) > 512)
       break;
    end    
    %}
    
    grad2 = Grad(x);
    y_k = grad2 - grad1;
    y_k = y_k';
    grad1 = grad2;
    
    % inputs for updating hessian matrix
    %{
    x1 = B * s_k';
    x2 = s_k * B;
    x3 = x2 * s_k';    
    x4 = y_k' * y_k;
    x5 = y_k' * s_k';
    x6 = x4/x5;
    x7 = (x1 * x2)/x3;
    % update Hessian matrix B^-1
    B = B - x7 + x6;
    H = inv(B);
    %}
     
    x1 = 1/(y_k'*s_k');
    if k ==1
        %H = (y_k'* s_k')/(y_k' * y_k) * eye(2);
    else
        H = (eye(2)- x1 * s_k * y_k) * H * (eye(2) - x1 * y_k' * s_k') ...
            + (x1 * (s_k * s_k'));
    end
    x_x(k) = x(1);
    x_y(k) = x(2);
    f_k(k) = f(x(1), x(2));    
    k = k+1;
    
end %while

if nargin > 5 && strcmp(plot,'plot')
    %plot graph of function and path
    
    
    rosenbrock_2d([x_0(1), x_0(2)],min(min(x_x, x_y)),max(max(x_x, x_y)));
    %test_function([X(1), X(2)],min(min(x_x, x_y)),max(max(x_x, x_y))) ;
    hold on
    plot3(x_x, x_y, f_k, 'r');
end
end
%end BFGS

%backtracking algorithm
function alpha = linesearch(a_k, beta, p_k, f, x)
 in = x + a_k * p_k';
 c1 = 1e-4;
 while feval(f, in(1), in(2)) > (feval(f, x(1), x(2)) + c1 * a_k * (p_k' * p_k))
    a_k = beta * a_k;
    in = x + a_k * p_k';
 end
 alpha = a_k;
end
% inexact line search - cubic interpolation method
function alpha = inexact(a_k, p_k, f, x, grad)
 in = x + a_k * p_k'; 
 c1 = 1e-4;
 n =1;
 f0(n) = feval(f, x(1), x(2));
 f1(n)= feval(f, in(1), in(2));
 
 alpha(n) = 1;
 %alpha(n) = 2*(f1(n) - f0(n))/(grad*p_k);
 
 
 while feval(f, in(1), in(2)) > (feval(f, x(1), x(2)) + c1 * a_k * (p_k' * p_k))
     n = n+1;
     if n == 2
         %quadratic interpolation find a_1
         alpha(n) = -a_k^2 * grad * p_k /(2*(f1(n-1)-f0(n-1)-grad * p_k*a_k));
     elseif n>2
         %cubic interpolation to find a_2 ...
         ab = [alpha(n-1)^3, alpha(n-1)^2; alpha(n-2)^3, alpha(n-2)^2];
         aa = inv(ab) * [f1(n-1)-alpha(n-1)*grad*p_k - f0(n-1); f1(n-2)-alpha(n-2)*grad*p_k - f0(n-2)];
         a = aa(1);
         b = aa(2);
         alpha(n) = (-b + sqrt(b^2 - 3*a*grad*p_k))/(3*a);
         
     end %if-else
    x = in;
    in = x + alpha(n) * p_k';
    f0(n) = f1(n-1);
    f1(n) = feval(f, in(1), in(2));
 
    if (abs((alpha(n)-alpha(n-1))/alpha(n-1))>0.9 || abs((alpha(n)-alpha(n-1))/alpha(n-1))<0.1)
     alpha(n) = 1/2 * alpha(n-1);
    end %if 
 
 end % while

 
 alpha = alpha(n);

end %inexact search

function out = Grad(x_k)
f_grad = @(x,y) [2*x - 400*x*(- x^2 + y) - 2, - 200*x^2 + 200*y];
%f_grad = @(x,y) [(x*sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/(2*abs(y - x + 47)^(1/2)) - sin(abs(y - x + 47)^(1/2)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(4*abs(x/2 + y + 47)^(1/2)), ...
%    - sin(abs(x/2 + y + 47)^(1/2)) - (x*sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/(2*abs(y - x + 47)^(1/2)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(2*abs(x/2 + y + 47)^(1/2))];
out = f_grad(x_k(1), x_k(2));
end

function out = Hess(x_k)
f_hess = @(x,y) [1200*x^2-400*y+2,-400*x;-400*x,200];
out = f_hess(x_k(1), x_k(2));
end