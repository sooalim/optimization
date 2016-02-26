f = @(x, y) -1*(y+47)*sin(sqrt(abs(x/2+(y+47))))-x*sin(sqrt(abs(x-(y+47)))); % eggholder function
k = 1;
result_CG = zeros(1, 6);
for i = 400:1:512
    for j =400:1:512
        try
            [a,f_k,x,y] = conjugate_gradient([i, j], f, 1e-6, 2000, 0.3, '');
        catch
            continue;
        end
        len = length(x);
        result_CG(k, :) = [i, j, x(len), y(len),f_k(len), len];
        k = k+1;
    end
end
T = array2table(result_CG);
T.Properties.VariableNames = {'i' 'j' 'x' 'y' 'f_k' 'N'};
rows = T.f_k == min(T.f_k);
minimum = T(rows, T.Properties.VariableNames)