% get data.
[XTrain, YTrain, anglesTrain] = digitTrain4DArrayData;
[XTest, YTest, anglesTest] = digitTest4DArrayData;


% get matrix dimensions and show sample image.
train_size = size(XTrain);
test_size = size(XTest);
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
net = patternnet(numel(class_names));
[net, tr] = train(net, XTrain, YTrain_OneHot);
% test the network.
predictions = net(XTest);