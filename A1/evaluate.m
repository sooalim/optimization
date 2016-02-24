function [ce, frac_correct] = evaluate(targets, y)
%    Compute evaluation metrics.
% to evaluate how good the classifier is
%    Inputs:
%        targets : N x 1 vector of binary targets. Values should be either 0 or 1.
%        y       : N x 1 vector of probabilities.
%    Outputs:
%        ce           : (scalar) Cross entropy. CE(p, q) = E_p[-log q]. Here we
%                       want to compute CE(targets, y).
%        frac_correct : (scalar) Fraction of inputs classified correctly.
% TODO: Finish this function

error(nargchk(2,2,nargin));
if (size(y) == size(targets))
    labels = y>=0.5;
    frac_correct = mean(targets==labels); % classification rate of y
    ce=sum(-targets.*log(y)-(1-targets).*log(1-y)); %cross entropy betewen target and output
end
