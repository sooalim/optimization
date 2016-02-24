function [X,F,Iters, X_i] = bfgs(N, X, gradToler, DxToler, MaxIter, myFx)
% Function bfgs performs multivariate optimization using the
% Broyden-Fletcher-Goldfarb-Shanno method.
%
% Input
%
% N - number of variables
% X - array of initial guesses
% gradToler - tolerance for the norm of the slopes
% DxToler - array of delta X tolerances
% MaxIter - maximum number of iterations
% myFx - name of the optimized function
%
% Output
%
% X - array of optimized variables
% F - function value at optimum
% Iters - number of iterations
%

B = eye(N,N)
bGoOn = true;
Iters = 0;
% calculate initial gradient
X_i(Iters+1, 1:5) = [Iters,X, helical(X)];
grad1 =  Grad(X, N, myFx);
grad1 = grad1';

while bGoOn

  Iters = Iters + 1;
  if Iters > MaxIter
    break;
  end
 
  P = -1 * B^-1 * grad1;
  P = P' / norm(P); % normalize direction vector P

  alpha = 1;
  alpha = linsearch(X, N, alpha, P, myFx);
  % calculate optimum X() with the given Lambda
  d = alpha * P;
  GRAD(Iters, 1:3) = grad1';
  P_i(Iters, 1:3) = P;
  X = X + d;
  X_i(Iters+1, 1:5) = [Iters,X, helical(X)];
  % get new gradient
  grad2 =  Grad(X, N, myFx);
  grad2 = grad2';
  g = grad2 - grad1;
  grad1 = grad2;

  % test for convergence
  for i = 1:N
    if abs(d(i)) > DxToler(i)
      break
    end
  end

  if norm(grad1) < gradToler
    break
  end

  d = d';
  x1 = g * g';
  x2 = g' * d;
  x3 = B * d;
  x4 = d' * B;
  x5 = d' * B * d;
  
  B = B + x1 / x2 - (x3 * x4) / x5
  
end

F = feval(myFx, X);
X_i(Iters+1, 1:5) = [Iters,X, helical(X)];

%plot graphs
subplot(2,1,1)
plot(X_i(:, 1), X_i(:, 2))
hold on
plot(X_i(:, 1), X_i(:, 3))
plot(X_i(:, 1), X_i(:, 4))
title('(BFGS: Quasi-Newton) Iteration number(k) vs x_i for i = 1, 2, 3')
legend('x1', 'x2', 'x3')
hold off
subplot(2,1,2)
plot(X_i(:, 1), X_i(:, 5))
title('Iteration number(k) vs f(x)')
legend('f(x)')
% end

function y = myFxEx(X, DeltaX, alpha, myFx)

  X = X + alpha * DeltaX;
  y = feval(myFx, X);
 
% end

function FirstDerivX = Grad(X, N, myFx)
%finite forward method for differentiation

for iVar=1:N
  xt = X(iVar);
  h = 0.01 * (1 + abs(xt));
  X(iVar) = xt + h;
  fp = feval(myFx, X);
  X(iVar) = xt - h;
  fm = feval(myFx, X);
  X(iVar) = xt;
  FirstDerivX(iVar) = (fp - fm) / 2 / h;
end

% end

function alpha = linsearch(X, N, alpha, D, myFx)

  MaxIt = 100;
  Toler = 0.1e-6;

  iter = 0;
  bGoOn = true;
  while bGoOn
    iter = iter + 1;
    if iter > MaxIt
      alpha = 0;
      break
    end

    h = 0.01 * (1 + abs(alpha));
    f0 = myFxEx(X, D, alpha, myFx);
    fp = myFxEx(X, D, alpha+h, myFx);
    fm = myFxEx(X, D, alpha-h, myFx);
    deriv1 = (fp - fm) / 2 / h;
    deriv2 = (fp - 2 * f0 + fm) / h ^ 2;
    if deriv2 == 0
      break
    end
    diff = deriv1 / deriv2;
    alpha = alpha - diff;
    if abs(diff) < Toler
      bGoOn = false;
    end
  end

% end
