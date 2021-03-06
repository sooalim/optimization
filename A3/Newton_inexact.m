function [x,f_k, x_x, x_y, step_length, z] = Newton_inexact(x_0, f, tol, MaxIter, beta, plot)

format long

if (nargin < 4) 
    tol = 1e-6;
end %if

k = 1; %initialize count
x= x_0;

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
    a_k = inexact(a_k, p_k, f, x, grad1);
    
    s_k = a_k * p_k; % x_(k+1) - x_k
    s_k = s_k'; 
    
    old_x = x;
    x = x + s_k;
    if(x(1) < -512 || x(2) < -512 || x(1) > 512 || x(2) > 512)
       break;
    end    
    step_length(k, 1:2) = [k, norm(x-old_x, 2)];
    
    B = Hess(x);
    H = inv(B);
    
    %storing iterates
    x_x(k) = x(1);
    x_y(k) = x(2);
    z(k, 1:3) = [k,x];
    f_k(k) = f(x(1), x(2));
    norm_grad(k) = norm(grad1, 2);
    
    k = k+1;
    
end %while

% plot graph of path and function if 'plot'
if nargin > 5 && strcmp(plot,'plot')
    %plot graph of function and path      
    %rosenbrock_2d([x_0(1), x_0(2)],min(min(x_x, x_y)),max(max(x_x, x_y)));
    test_function([x_0(1), x_0(2)],min(min(x_x, x_y)),max(max(x_x, x_y))) ;
    hold on
    plot3(x_x, x_y, f_k, 'r');
    scatter3(x_x, x_y, f_k, 'b*');
end
end
%end Newton

% inexact line search - quadratic & cubic interpolation method
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
         aa = ab^-1 * [f1(n-1)-alpha(n-1)*grad*p_k - f0(n-1); f1(n-2)-alpha(n-2)*grad*p_k - f0(n-2)];
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

%hard coded gradient and hessian function 
function out = Grad(x_k)
%f_grad = @(x,y) [2*x - 400*x*(- x^2 + y) - 2, - 200*x^2 + 200*y];
f_grad = @(x,y) [(x*sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/(2*abs(y - x + 47)^(1/2)) - sin(abs(y - x + 47)^(1/2)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(4*abs(x/2 + y + 47)^(1/2)), ...
    - sin(abs(x/2 + y + 47)^(1/2)) - (x*sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/(2*abs(y - x + 47)^(1/2)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(2*abs(x/2 + y + 47)^(1/2))];
out = f_grad(x_k(1), x_k(2));
end

function out = Hess(x_k)
%rosenbrock hessian
%f_hess = @(x,y) [1200*x^2-400*y+2,-400*x;-400*x,200];

%egg_holder hessian

f_hess = @(x,y) [(sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/abs(y - x + 47)^(1/2) + (sign(x/2 + y + 47)^2*sin(abs(x/2 + y + 47)^(1/2))*(y + 47))/(16*abs(x/2 + y + 47)) - (x*dirac(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/abs(y - x + 47)^(1/2) - (dirac(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(4*abs(x/2 + y + 47)^(1/2)) + (x*sign(y - x + 47)^2*cos(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)^(3/2)) + (x*sign(y - x + 47)^2*sin(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)) + (sign(x/2 + y + 47)^2*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(16*abs(x/2 + y + 47)^(3/2)), ...
(sign(x/2 + y + 47)^2*sin(abs(x/2 + y + 47)^(1/2))*(y + 47))/(8*abs(x/2 + y + 47)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2)))/(4*abs(x/2 + y + 47)^(1/2)) - (sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/(2*abs(y - x + 47)^(1/2)) + (x*dirac(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/abs(y - x + 47)^(1/2) - (dirac(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(2*abs(x/2 + y + 47)^(1/2)) - (x*sign(y - x + 47)^2*cos(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)^(3/2)) - (x*sign(y - x + 47)^2*sin(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)) + (sign(x/2 + y + 47)^2*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(8*abs(x/2 + y + 47)^(3/2)); ...
(sign(x/2 + y + 47)^2*sin(abs(x/2 + y + 47)^(1/2))*(y + 47))/(8*abs(x/2 + y + 47)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2)))/(4*abs(x/2 + y + 47)^(1/2)) - (sign(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/(2*abs(y - x + 47)^(1/2)) + (x*dirac(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/abs(y - x + 47)^(1/2) - (dirac(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(2*abs(x/2 + y + 47)^(1/2)) - (x*sign(y - x + 47)^2*cos(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)^(3/2)) - (x*sign(y - x + 47)^2*sin(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)) + (sign(x/2 + y + 47)^2*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(8*abs(x/2 + y + 47)^(3/2)), ... 
(sign(x/2 + y + 47)^2*sin(abs(x/2 + y + 47)^(1/2))*(y + 47))/(4*abs(x/2 + y + 47)) - (sign(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2)))/abs(x/2 + y + 47)^(1/2) - (x*dirac(y - x + 47)*cos(abs(y - x + 47)^(1/2)))/abs(y - x + 47)^(1/2) - (dirac(x/2 + y + 47)*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/abs(x/2 + y + 47)^(1/2) + (x*sign(y - x + 47)^2*cos(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)^(3/2)) + (x*sign(y - x + 47)^2*sin(abs(y - x + 47)^(1/2)))/(4*abs(y - x + 47)) + (sign(x/2 + y + 47)^2*cos(abs(x/2 + y + 47)^(1/2))*(y + 47))/(4*abs(x/2 + y + 47)^(3/2))];
 

out = f_hess(x_k(1), x_k(2));
end