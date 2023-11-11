% refresh session.
clf;

% load the data files.
dataColumnNames = {'Packing' 'VoroMQA' 'SBROD' '3DCNN' 'ProQ2' 'ProQ3' 'ProteinPro' 'GDT', 'predGDT'};


% train data.
train_file = 'data/decoy_train_features.txt';
train_data = readmatrix(train_file, 'NumHeaderLines', 1);
XTrain = train_data(:, 2:8)';
TTrain = train_data(:, 9)';


% test data.
test_file = 'data/decoy_test_features.txt';
test_data = readmatrix(test_file, 'NumHeaderLines', 1);
XTest = test_data(:, 2:8)';
TTest = test_data(:, 9)';


% NN 1. -----
% train network.
net = feedforwardnet(10);
[net, tr] = train(net, XTrain, TTrain);
figure(1);
plotperform(tr);

% test network.
YTrain = net(XTrain);
YTest_nn = net(XTest);
test_mae = mae(YTest_nn, TTest);
fprintf("NN Test MAE: %f\n", test_mae);

% plot error histogram.
error = TTest - YTest_nn;
figure('Name', 'NN 1');
ploterrhist(error);
annotation('textbox', [0.2 0.5 0.3 0.3], 'String', sprintf("NN 1 MAE = %f", test_mae)', 'FitBoxToText', 'on');


% NN 2. -----
% train network.
net = feedforwardnet(20);
[net, tr] = train(net, XTrain, TTrain);
figure(1);
plotperform(tr);

% test network.
YTrain = net(XTrain);
YTest_nn = net(XTest);
test_mae = mae(YTest_nn, TTest);
fprintf("NN Test MAE: %f\n", test_mae);

% plot error histogram.
error = TTest - YTest_nn;
figure('Name', 'NN 2');
ploterrhist(error);
annotation('textbox', [0.2 0.5 0.3 0.3], 'String', sprintf("NN 2 MAE = %f", test_mae)', 'FitBoxToText', 'on');


% NN 3. -----
% train network.
net = feedforwardnet([15 5]);
[net, tr] = train(net, XTrain, TTrain);
figure(1);
plotperform(tr);

% test network.
YTrain = net(XTrain);
YTest_nn = net(XTest);
test_mae = mae(YTest_nn, TTest);
fprintf("NN Test MAE: %f\n", test_mae);

% plot error histogram.
error = TTest - YTest_nn;
figure('Name', 'NN 3');
ploterrhist(error);
annotation('textbox', [0.2 0.5 0.3 0.3], 'String', sprintf("NN 3 MAE = %f", test_mae)', 'FitBoxToText', 'on');


% NN 4. -----
% train network.
net = feedforwardnet([20 10]);
[net, tr] = train(net, XTrain, TTrain);
figure(1);
plotperform(tr);

% test network.
YTrain = net(XTrain);
YTest_nn = net(XTest);
test_mae = mae(YTest_nn, TTest);
fprintf("NN Test MAE: %f\n", test_mae);

% plot error histogram.
error = TTest - YTest_nn;
figure('Name', 'NN 4');
ploterrhist(error);
annotation('textbox', [0.2 0.5 0.3 0.3], 'String', sprintf("NN 4 MAE = %f", test_mae)', 'FitBoxToText', 'on');


% plot correlations.
figure(8);
corrplot(XTrain', 'varNames', dataColumnNames(1:7));

figure(9);
corrplot([XTrain; TTrain; YTrain]', 'varNames', dataColumnNames(1:9));


% try and fit a linear model.

% One feature, linear model.
linear_mdl = fitlm(XTrain(1, :)', TTrain');
YTest_linear = predict(linear_mdl, XTest(1, :)');
test_mae_linear = mae(YTest_linear', TTest);
fprintf("Linear Model Test MAE: %f\n", test_mae_linear);

% Six feature, linear model.
% linear_mdl = fitlm(XTrain(1:6, :)', TTrain');
% YTest_linear = predict(linear_mdl, XTest(1:6, :)');
% test_mae_linear = mae(YTest_linear', TTest);
% fprintf("Linear Model Test MAE: %f\n", test_mae_linear);

% plot error histogram for linear model.
figure(10);
plotResiduals(linear_mdl);
annotation('textbox', [0.2 0.5 0.3 0.3], 'String', sprintf("MAE = %f", test_mae_linear)', 'FitBoxToText', 'on');