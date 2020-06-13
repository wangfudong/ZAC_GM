%addpath(genpath(pwd))
model = textread('E:\computer_vision\code\point-reg\gmmreg-master\data\fish_data\fish_X.txt');
scene = textread('E:\computer_vision\code\point-reg\gmmreg-master\data\fish_data\fish_Y.txt');
[config] = initialize_config(model,scene,'tps');
[param, transformed_model, history, config] = gmmreg_L2(config);
figure(1)
DisplayPoints2D(model,scene);
figure(2)
DisplayPoints2D(transformed_model,scene);