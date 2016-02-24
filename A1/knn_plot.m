%function knn_plot(ti, tt, vi, vt) 
% knn_plot: Script that runs kNN for different values of k in {1, 3, 5, 7, 9} 
%           and plots the classification rate on the validation set as a
%           function of k.
clear all;
close all;
load mnist_train;
load mnist_test;
load mnist_valid;


k=1:2:9;
ce = ones(size(k,2),1);

for i=1:size(k,2)
    t=run_knn(k(i), train_inputs, train_targets, valid_inputs);
    ce(i) = mean(t==valid_targets);
    %classification rate for validation set
    
    d=run_knn(k(i), train_inputs, train_targets, test_inputs);
    ce_test(i) = mean(d==test_targets);
    %classification rate for test set
end

plot(k,ce, '-o');
hold;
grid on;
plot(k,ce_test, '-o');
legend('validation set','test set');
title('Classification rate for different k');
ylabel('Classification rate');
xlabel('k');
axis([0,max(k)+1,0.8,1.0]);
hold;





