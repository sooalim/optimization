function [f, df, y] = logistic(weights, data, targets, hyperparameters)
% Calculate log likelihood and derivatives with respect to weights.
%
% Note: N is the number of examples and 
%       M is the number of features per example.
%
% Inputs:
% 	weights:    (M+1) x 1 vector of weights, where the last element
%               corresponds to bias (intercepts).
% 	data:       N x M data matrix where each row corresponds 
%               to one data point.
%	targets:    N x 1 vector of binary targets. Values should be either 0 or 1.
%   hyperparameters: The hyperparameter structure
%
% Outputs:
%	f:             The scalar error value.
%	df:            (M+1) x 1 vector of derivatives of error w.r.t. weights.
%   y:             N x 1 vector of probabilities. This is the output of the classifier.
%
%TODO: finish this function
% compute E(w,b) and dE(w,b)/dw
N = size(data, 1) %number of training data
newX = [data';ones(1,N)]';
wx = weights'*newX'; %compute w_i*x_i+b 
                                  %exponential of logistic function
                                  %append original data with ones to
                                  %accomodate for b, the bias term (= w_0)
%D = size(data, 2) %number of featuresof training data
f = 0;
for i=1:N;    
    s = sigmoid(wx(i));
    dlikelihood = -targets.*wx(i)+wx(i)-log(1/s);
    f = f + dlikelihood;
    df = -data(i, :)* targets + s;
end

N = length(targets); % number of examples
f = evaluate(targets, y);
df = 
