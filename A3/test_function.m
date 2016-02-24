function varargout = test_function(X, r_x, r_y)
% Rosenbrock function
%
%   eggholder([x1, x2]) returns the value of the value of the 
%   two dimensional eggholder function at the specified points. All [xi] may be 
%   vectors. 
%
%   The global minimum is
%
%                   f(x1, x2) = ...


    % if no input is given, return dimensions, bounds and minimum
    f = @(x, y) -1*(y+47)*sin(sqrt(abs(x/2+(y+47))))-x*sin(sqrt(abs(x-(y+47)))); % eggholder function
    
    if (nargin == 0)
        varargout{1} = [512, 404.2319];  % #solution
        varargout{2} = f(512, 404.2319);
        
    end
    if(nargin > 0)
        
        % split input vector X into x1, x2, x3
        if size(X, 1) == 2
            x1 = X(1, :);  x2 = X(2, :);
        else
            x1 = X(:, 1);  x2 = X(:, 2);
        end
        
        % output function value
        varargout{1} = f(x1, x2);
    end
     
    if (nargin == 3)
       syms x y
       ezsurfc(f, [r_x, r_y, r_x, r_y], 100)
       xlabel('x'); ylabel('y'); zlabel('function'); colorbar 
       text(512, 404.2319 ,1000, '\downarrow function minimal point(512, 404.2319)');
       colormap(jet)
       colorbar
    end

end