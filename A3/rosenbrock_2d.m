function varargout = rosenbrock_2d(X, r_x, r_y)
% Rosenbrock function
%
%   rosenbrock([x1, x2]) returns the value of the value of the 
%   two dimensional rosenbrock function at the specified points. All [xi] may be 
%   vectors. 
%
%   The global minimum is
%
%                   f(x1, x2) = f(1, 0) = 
%
% Please report bugs and inquiries to: 

    % if no input is given, return dimensions, bounds and minimum
    f = @(x, y) 100*(y - x.^2).^2 + (1-x).^2;
    if (nargin == 0)
        varargout{1} = [1, 0];  % # dims
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
       s = linspace(r_x,r_y,40);
       t = s;
       colormap(white);
       [x, y] = meshgrid(s, t);
       surfc(x, y, 100*(y - x.^2).^2 + (1-x).^2) 
       xlabel('x'); ylabel('y'); zlabel('function'); %colorbar
       t = text(1, 1 , f(1,1), '\leftarrow function minimal point');
       t.Color = 'red';
       t.FontSize = 15;
    end

end