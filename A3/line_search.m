function [alpha] = line_search(fcn,x,ds,finf,sigma,guess,maxit);

% LINE_SEARCH: A line search descent method under improved Wolfe-Powell conditions: 
% it imposes a more stringent two-sided test on the slope, which allows the user to restrict 
% the set of acceptable points to a small neighbourhood of a local minimizer
% by reducing 'sigma'. Moreover, cubic interpolation is applied. Based on:
%    Wolfe (1968): "Convergence conditions for ascent methods",
%    Powell (1976): "Some global convergence properties of a variable
%                    metric algorithm for minimization without exact line searches", and 
%    Fletcher (2006): "Practical Methods of Optimization".
%
% inputs:
%      required:
%       fcn: the name of a matlab function that must return the value
%       of the function to be minimized 'f' and its first-derivative 'df'
%       evaluated at a given point 'x', a given search direction 'ds',
%       and a given guess "alpha".
%       x: a given point of the domain of 'f'.
%       ds: a given search direction.
%       finf: a lower bound on f. Any f<finf would be optimal.
%            examples: firms' optimization problem: finf=0, 
%                      nonlinear least squares problem: finf=0.
%      optional:
%       sigma: 0<=sigma<=1, a sigma below or equal to 0.1 is equivalent
%             to perform an exact line search. Default: sigma=0.9.
%       guess: Default: guess=1.
%       maxit: maximum number of iterations. Default: maxit=15.
%
% SEE EXAMPLE: line_search_ex.txt
%
% Angel Luis Lopez. 
% email: angel.lopez@univ-tlse1.fr

      ro = 0.01;
      t1 = 9;
      t2 = 0.1;
      t3 = 0.5;
      
      flag = 0;
      
      alpha(1) = 0;
      
      if nargin<5, sigma = 0.9; end
      if nargin<6, alpha(2) = 1; else alpha(2) = guess; end
      if nargin<7, maxit = 15; end

      [f0,df0] = feval(fcn,x,ds,alpha(1));
      
      if df0>0 error('descent condition is not satisfied: df>0!'); end
      
      f(1)=f0; df(1)=df0;
      
      u = (finf-f0)/(ro*df0);
      
      % BRACKETING PHASE
      
      for i=2:maxit
          
          [f(i),df(i)] = feval(fcn,x,ds,alpha(i));
          
          if f(i) <= finf flag=1; break, end % terminate 'line search'
          
          if f(i) > f0 + alpha(i)*ro*df0 | f(i) >= f0
              
                a(i) = alpha(i-1); b(i) = alpha(i); break, end % terminate 'bracketing phase'
               
          if abs(df(i)) <= -sigma*df0' flag=1; break % terminate 'line search'
              
          elseif df(i) >= 0
              
                a(i) = alpha(i); b(i) = alpha(i-1); break, end % terminate 'bracketing phase'
            
          if u <= 2*alpha(i) - alpha(i-1)
              
              alpha(i+1) = u;
              
          else
              
              % interval
              
              intval = [2*alpha(i)-alpha(i-1), min(u,alpha(i)+t1*(alpha(i)-alpha(i-1)))];
              
              %%% mapping on to [0,1]
              
              map = [alpha(i-1) alpha(i)];
              dfm(i-1) = (map(2)-map(1))*df(i-1);
              dfm(i) = (map(2)-map(1))*df(i);
              
              % parameters of the Hermite interpolating cubic
              
              c1 = f(i-1);
              c2 = dfm(i-1);
              c3 = 3*(f(i)-f(i-1)) - 2*dfm(i-1) - dfm(i);
              c4 = dfm(i-1) + dfm(i) - 2*(f(i)-f(i-1));
              
              % interval: alpha = a + z(b-a); alpha is in intval
              
              zmin = (intval(1)-map(1))/(map(2)-map(1));
              zmax = (intval(2)-map(1))/(map(2)-map(1));
              
              % min c(z) = c1+c2z+c3z^2+c4z^3,
              % where z belongs to [zmin,zmax]
              
              mincubic = (-c3+sqrt(c3^2-3*c4*dfm(i-1)))/(3*c4);
              
              if mincubic<zmin, mincubic=zmin;
              elseif mincubic>zmax, mincubic=zmax; end
       
              z=[zmin;zmax;mincubic];
              
              c = c1+c2*z+c3*z.^2+c4.*z.^3;
              
              [cmin,indc] = min(c);
              
              z = z(indc); % optimal z
              
              % New value for alpha
              
              alpha(i+1) = map(1) + z*(map(2)-map(1));
              
          end
          
      end
      
      if flag alpha = alpha (i); return
      
     % SECTIONING PHASE 
     
      else 
          
          for k=i:maxit
              
              % mapping [a,b] on to [0,1]
                
                [fa(k),dfa(k)] = feval(fcn,x,ds,a(k));

                [fb(k),dfb(k)] = feval(fcn,x,ds,b(k));
                
                intval = [a(k)+t2*(b(k)-a(k)),b(k)-t3*(b(k)-a(k))];
                
                %%% mapping on to [0,1]
              
                map = [a(k) b(k)];
                
                dfam(k) = (map(2)-map(1))*dfa(k);
                dfbm(k) = (map(2)-map(1))*dfb(k);
              
              % parameters of the Hermite interpolating cubic
              
                c1 = fa(k);
                c2 = dfam(k);
                c3 = 3*(fb(k)-fa(k)) - 2*dfam(k) - dfbm(k);
                c4 = dfam(k) + dfbm(k) - 2*(fb(k)-fa(k));
                
              % interval: alpha = a + z(b-a); alpha is in intval
              
                zmin = (intval(1)-map(1))/(map(2)-map(1));
                zmax = (intval(2)-map(1))/(map(2)-map(1));
                
              % min c(z) = c1+c2z+c3z^2+c4z^3,
              % where z belongs to [zmin,zmax]
              
                mincubic = (-c3+sqrt(c3^2-3*c4*dfam(k)))/(3*c4); 
               
              
                if mincubic<zmin, mincubic=zmin;
                elseif mincubic>zmax, mincubic=zmax; end
       
                z=[zmin;zmax;mincubic];
                
                
                c = c1+c2*z+c3*z.^2+c4.*z.^3;
              
                [cmin,indc] = min(c);
              
                z = z(indc); % ---> optimal z             
                
                % New value for alpha
              
                alpha(k) = map(1) + z*(map(2)-map(1));
                
                [f(k),df(k)] = feval(fcn,x,ds,alpha(k));
                
                if f(k) > f0+ro*alpha(k)*df0 | f(k) >= fa(k)
                    
                    a(k+1) = a(k); b(k+1) = alpha(k);
                    
                else
                    
                    if abs(df(k)) <= -sigma*df0 alpha=alpha(k); return, end; % terminate 'line search'
                    
                    a(k+1) = alpha(k);
                    
                    if (b(k)-a(k))*df(k) >= 0
                        
                        b(k+1)=a(k);
                        
                    else 
                        
                        b(k+1)=b(k);
                        
                    end
                    
                    if (a(k)-alpha(k))*dfa(k) <= sqrt(eps)
                        
                        warning('Potential Round-off error, no progress is posible in the line search'); break
                        
                    end
                    
                end
             
            end
            
      end
                
              