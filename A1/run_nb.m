% Learn a Naive Bayes classifier on the digit dataset, evaluate its
% performance on training and test sets, then visualize the mean and variance
% for each class.

load mnist_train;
load mnist_test;

% Add your code here (it should be less than 10 lines)

[log_prior, class_mean, class_var] = train_nb(train_inputs, train_targets); %train model

plot_digits(class_mean) %visualization of mean vector for both classes
plot_digits(class_var) %visualization of variance vector for both classes

%compute accuracy on test data
[test_prediction, test_accuracy] = test_nb(test_inputs, test_targets, log_prior, class_mean, class_var); 

%compute accuracy on train data
[tr_prediction, tr_accuracy] = test_nb(train_inputs, train_targets, log_prior, class_mean, class_var); 

disp(['test_accuracy: ' num2str(test_accuracy)]) ;
disp(['train_accuracy: ' num2str(tr_accuracy)]) ;