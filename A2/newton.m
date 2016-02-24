function [X,F,Iters, X_in] = newton(N, X, gradToler, XToler, MaxIter, myFx)
% Function NEWTON_OPT performs multivariate optimization using the
% Newton's method.
%
% Input
%
% N - number of variables
% X - array of initial guesses
% gradToler - tolerance for the norm of the slopes
% MaxIter - maximum number of iterations
% myFx - name of the optimized function
%
% Output
%
% X - array of optimized variables
% F - function value at optimum
% Iters - number of iterations
%

bGoOn = true;
Iters = 0;
X_in(Iters+1, 1:5) = [Iters, X(1), X(2), X(3),feval(myFx, X)];
while bGoOn

  Iters = Iters + 1
  if Iters > MaxIter
    break;
  end

  g = FirstDerivatives(X, N, myFx);
  fnorm = norm(g);
  if fnorm < gradToler
    break;
  end
  J = SecondDerivatives(X, N, myFx);
  DeltaX = g / J;

  X = X - DeltaX
  bStop = true;
  for i=1:N
    if abs(DeltaX(i)) > XToler(i)
      bStop = false;
    end
  end

  bGoOn = ~bStop;
  fx = feval(myFx, X);
  X_in(Iters+1, 1:5) = [Iters, X(1), X(2), X(3), fx]
end

F = feval(myFx, X);
X_in(Iters+1, 1:5) = [Iters, X(1), X(2), X(3), F]

%plot graphs
subplot(2,1,1)
plot(X_in(:, 1), X_in(:, 2))
hold on
plot(X_in(:, 1), X_in(:, 3))
plot(X_in(:, 1), X_in(:, 4))
title('(Newton direction) Iteration number(k) vs x_i for i = 1, 2, 3')
legend('x1', 'x2', 'x3')
hold off
subplot(2,1,2)
plot(X_in(:, 1), X_in(:, 5))
title('Iteration number(k) vs f(x)')
legend('f(x)')
% end

function FirstDerivX = FirstDerivatives(X, N, myFx)

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

function SecondDerivX = SecondDerivatives(X, N, myFx)

for i=1:N
  for j=1:N
    % calculate second derivative?
    if i == j
      f0 = feval(myFx, X);
      xt = X(i);
      hx = 0.01 * (1 + abs(xt));
      X(i) = xt + hx;
      fp = feval(myFx, X);
      X(i) = xt - hx;
      fm = feval(myFx, X);
      X(i) = xt;
      y = (fp - 2 * f0 + fm) / hx ^ 2;
    else
      xt = X(i);
      yt = X(j);
      hx = 0.01 * (1 + abs(xt));
      hy = 0.01 * (1 + abs(yt));
      % calculate fpp;
      X(i) = xt + hx;
      X(j) = yt + hy;
      fpp = feval(myFx, X);
      % calculate fmm;
      X(i) = xt - hx;
      X(j) = yt - hy;
      fmm = feval(myFx, X);
      % calculate fpm;
      X(i) = xt + hx;
      X(j) = yt - hy;
      fpm = feval(myFx, X);
      % calculate fmp
      X(i) = xt - hx;
      X(j) = yt + hy;
      fmp = feval(myFx, X);
      X(i) = xt;
      X(j) = yt;
      y = (fpp - fmp - fpm + fmm) / (4 * hx * hy);
    end
    SecondDerivX(i,j) = y;
  end
end
%end
