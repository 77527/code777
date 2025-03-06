clear;
close all;
clc;

NodeAmountArr=[100 150 200 250 300 350 400] ; 
BorderLength = 100; 
times = 1; 

for run=1:times
for index = 1:length(NodeAmountArr)
NodeAmount=NodeAmountArr(index);
%Nodes were randomly generated within the regional range
AreaC = BorderLength.*rand(2,NodeAmount);%[x1,...,xn;y1,...,yn;];
BeaconAmount = NodeAmount * 0.2; 
UnAmount = NodeAmount - BeaconAmount; 
R = 30; 
%Add ordnumbers to each point
data = [(1:NodeAmount);AreaC];
%Beacon coordinate information
BeaconData = data(2:3,1:BeaconAmount);
UnKnownData = data(2:3,BeaconAmount+1:end);

%% algorithm1£ºoriginal Dvhop
[X,d]=Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData);

%% algorithm2£ºDvhop for multiple communication radii and weighted hop distances
[X_w,d_w]=MRW_Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData);

%% algorithm3£ºSSAPSO-Dvhop
%% Set the parameters
XSSA = [];
for i = 1:UnAmount
    dim = 2; 
    pop = 20;
    MaxIter = 20;
    ub = BorderLength; 
    lb = 0;
    fobj = @(x) fun(x,UnAmount,BeaconAmount,d,BeaconData,i); 
    [Best_score,Best_pos,SSA_curve(i,:)]=SSAPSO(pop,MaxIter,lb,ub,dim,fobj); 
    XSSA = [XSSA;Best_pos];
end
XSSA =XSSA';

%% algorithm4£ºMRW-SSAPSO-Dvhop
%% Set the parameters
XWSSA = [];
for i = 1:UnAmount
    dim = 2; 
    pop = 20;
    MaxIter = 20;
    ub = BorderLength; 
    lb = 0;%
    fobj = @(x) fun(x,UnAmount,BeaconAmount,d_w,BeaconData,i);
    [Best_score_WSSA,Best_pos_WSSA,WSSA_curve(i,:)]=SSAPSO(pop,MaxIter,lb,ub,dim,fobj); 
    XWSSA = [XWSSA;Best_pos_WSSA];
end
XWSSA =XWSSA';

%Calculate the positioning error
clear error errorSSA error_WSSA
for i=1:UnAmount
    error(1,i)=((X(1,i)-UnKnownData(1,i))^2+(X(2,i)-UnKnownData(2,i))^2)^0.5;
    errorSSA(1,i) = ((XSSA(1,i)-UnKnownData(1,i))^2+(XSSA(2,i)-UnKnownData(2,i))^2)^0.5;
    error_WSSA(1,i) = ((XWSSA(1,i)-UnKnownData(1,i))^2+(XWSSA(2,i)-UnKnownData(2,i))^2)^0.5;
end


Accuracy(run,index)=sum(error)/(UnAmount*R);
AccuracySSA(run,index)=sum(errorSSA)/(UnAmount*R);
AccuracyWSSA(run,index)=sum(error_WSSA)/(UnAmount*R);
end
run
end

mean_Accuracy=sum(Accuracy,1)/times;
mean_AccuracySSA=sum(AccuracySSA,1)/times;
mean_AccuracyWSSA=sum(AccuracyWSSA,1)/times;

figure(1)
plot(NodeAmountArr(1:index),mean_Accuracy(1,1:index),'k-p','MarkerFaceColor','y')
hold on
plot(NodeAmountArr(1:index),mean_AccuracySSA(1,1:index),'k-^','MarkerFaceColor','g')
% hold on
% plot(NodeAmountArr(1:index),mean_AccuracyWSSA(1,1:index),'k-o','MarkerFaceColor','r')
xlabel('total number of nodes');
ylabel('average positioning error');
title('positioning error comparison(different total number of nodes)')
legend('Dvhop','SSAPSO-Dvhop')