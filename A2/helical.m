function fx = helical(x)
  if ( 0.0 < x(1) )
    theta = atan( x(2)/ x(1) ) / 2.0 / pi;
  elseif ( x(1) == 0.0 )
    theta = 0.25;
  elseif ( x(1) < 0.0 )
    theta = 0.5 + atan( x(2) / x(1) ) / 2.0 / pi;
  end
  fx1 = 10*(x(3) - 10.0 * theta);
  fx2 = 10*(sqrt ( x(1) * x(1) + x(2) * x(2) ) - 1);
  fx3 = x(3);

  fx =  fx1 * fx1 ...
     +  fx2 * fx2 ...
     +  fx3 * fx3;

  return
end

