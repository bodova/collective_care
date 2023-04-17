function[Split,idtemp] = Inference_TrainTestSplits(Nsplits,Ptrain);

Split = [];                                         % cell objects

idx = 1;
idtemp = 1;
[repTrain,repTest] = getrep_traintest(Ptrain);
Split{1,idx} = repTrain; 
Split{2,idx} = repTest;
