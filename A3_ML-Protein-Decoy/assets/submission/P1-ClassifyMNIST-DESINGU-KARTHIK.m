clf;

% get data.
[XTrain, YTrain, anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest] = digitTest4DArrayData;


% get matrix dimensions and show sample image.
train_size = size(XTrain);
test_size = size(XTest);
figure(1);
imshow(XTrain(:, :, :, 1));


% flatten the train and test arrays.
XTrain = reshape(XTrain, prod(train_size(1:3), "all"), train_size(4));
XTest = reshape(XTest, prod(test_size(1:3), "all"), test_size(4));
% update train and test sizes.
train_size = size(XTrain);
test_size = size(XTest);


% one-hot encode the Y vectors.
[class_ids, class_names] = grp2idx(YTrain);
YTrain_OneHot = zeros(numel(class_names), train_size(2));
for i = 1:train_size(2)
	YTrain_OneHot(class_ids(i), i) = 1;
end

YTest_OneHot = zeros(numel(class_names), test_size(2));
for i = 1:test_size(2)
	YTest_OneHot(class_ids(i), i) = 1;
end


% train a neural network; test it.
net = patternnet(10);
[net, tr] = train(net, XTrain, YTrain_OneHot);

disp(net.layers{1}.transferFcn);

% test the network; obtain predictions.
predictions = net(XTest)';
% compute test accuracy.
[max_values, predicted_num] = max(predictions, [], 2);
YTest_labels = double(YTest);
test_accuracy = sum(predicted_num == YTest_labels)./numel(YTest)*100;
fprintf("Test Accuracy: %f\n", test_accuracy);

% plot a confusion matrix on test predictions.
% figure(2);
% plotconfusion(YTest_labels, predicted_num);
figure(2);
confusionchart(YTest_labels, predicted_num);


% other neural net configurations.
activation_functions_c = {'compet', 'softmax', 'logsig', 'tansig', 'netinv'};
hidden_layers_c = {[25], [50], [25 15], [50 25], [50 25 15]};
accuracies = zeros(length(activation_functions_c), length(hidden_layers_c));

for j = 1:length(activation_functions_c)
	fprintf("Training for %s ...\n", activation_functions_c{j});

	for i = 1:length(hidden_layers_c)

		% init network.
		hidden_layers = hidden_layers_c{i};
		net = patternnet(hidden_layers);
		net.trainParam.showWindow = 0;  % don't show training window.
		[net, tr] = train(net, XTrain, YTrain_OneHot);

		% set activation function for each layer.
		for k = 1:length(hidden_layers)
			net.layers{k}.transferFcn = activation_functions_c{j};
		end

		% save results.
		predictions = net(XTest)';
		[max_values, predicted_num] = max(predictions, [], 2);
		accuracies(j, i) = sum(predicted_num == YTest_labels)./numel(YTest)*100;
	end
end
disp(accuracies);