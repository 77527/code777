clear;
close all;
clc;

radiusArr=[20 24 28 32 36 40]; % 通信半径数组
BorderLength = 100; %区域边界范围
NodeAmount = 250; %总的节点数
%在区域范围内随机生成节点，即总节点数NodeAmount个坐标
AreaC = BorderLength.*rand(2,NodeAmount);%[x1,...,xn;y1,...,yn;];
times = 1; %运行次数

for run=1:times
for index = 1:length(radiusArr)
BeaconAmount = 20; %锚节点数
UnAmount = NodeAmount - BeaconAmount; %未知节点数
R = radiusArr(index); %通信距离
%为每个点添加序号，如第1，2，3。放在第1行
data = [(1:NodeAmount);AreaC];
%信标坐标信息
BeaconData = data(2:3,1:BeaconAmount);%提取2，3行存放的坐标
UnKnownData = data(2:3,BeaconAmount+1:end);%提取剩下的坐标为未知节点坐标

%% 算法1：原始Dvhop
[X,d]=Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData);

%% 算法2：多通信半径和跳距加权的Dvhop
[X_w,d_w]=MRW_Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData);

%% 算法3：SSAPSO-Dvhop
%% 设置松鼠算法参数
XSSA = [];
for i = 1:UnAmount
    dim = 2; %dim是变量的数量（问题的维度）
    pop = 40;%种群数量
    MaxIter = 50;%最大迭代次数
    ub = BorderLength; %坐标最大范围
    lb = 0;%坐标最小范围
    fobj = @(x) fun(x,UnAmount,BeaconAmount,d,BeaconData,i);%适应度函数设置 
    [Best_score,Best_pos,SSA_curve(i,:)]=SSAPSO(pop,MaxIter,lb,ub,dim,fobj); %开始优化
    XSSA = [XSSA;Best_pos];
end
XSSA =XSSA';

%% 算法4：MRW-SSAPSO-Dvhop
%% 设置松鼠优化算法参数
XWSSA = [];
for i = 1:UnAmount
    dim = 2; %dim是变量的数量（问题的维度）
    pop = 20;%种群数量
    MaxIter = 20;%最大迭代次数
    ub = BorderLength; %坐标最大范围
    lb = 0;%坐标最小范围
    fobj = @(x) fun(x,UnAmount,BeaconAmount,d_w,BeaconData,i);%适应度函数设置，用的是d_w
    [Best_score_WSSA,Best_pos_WSSA,WSSA_curve(i,:)]=SSAPSO(pop,MaxIter,lb,ub,dim,fobj); %开始优化
    XWSSA = [XWSSA;Best_pos_WSSA];
end
XWSSA =XWSSA';

%计算误差
clear error errorSSA error_w errorWSSA
for i=1:UnAmount
    error(1,i)=((X(1,i)-UnKnownData(1,i))^2+(X(2,i)-UnKnownData(2,i))^2)^0.5;
    errorSSA(1,i) = ((XSSA(1,i)-UnKnownData(1,i))^2+(XSSA(2,i)-UnKnownData(2,i))^2)^0.5;
    errorWSSA(1,i) = ((XWSSA(1,i)-UnKnownData(1,i))^2+(XWSSA(2,i)-UnKnownData(2,i))^2)^0.5;
end

%基础Dvhop的误差计算
Accuracy(run,index)=sum(error)/(UnAmount*R);%归一化定位误差
%SSA_Dvhop的误差计算
AccuracySSA(run,index)=sum(errorSSA)/(UnAmount*R);
%MRW-SSA-hop的误差计算
AccuracyWSSA(run,index)=sum(errorWSSA)/(UnAmount*R);
end
run
end

mean_Accuracy=sum(Accuracy,1)/times;
mean_AccuracySSA=sum(AccuracySSA,1)/times;
mean_AccuracyWSSA=sum(AccuracyWSSA,1)/times;

figure(1)
plot(radiusArr(1:index),mean_Accuracy(1,1:index),'k-p','MarkerFaceColor','y')
hold on
plot(radiusArr(1:index),mean_AccuracySSA(1,1:index),'k-^','MarkerFaceColor','b')
% hold on
% plot(radiusArr(1:index),mean_AccuracyWSSA(1,1:index),'k-o','MarkerFaceColor','r')
xlabel('communication radius');
ylabel('average positioning error');
title('positioning error comparison(different communication radius)')
legend('Dvhop','SSAPSO-Dvhop')
% legend('原始Dvhop','SSA-Dvhop','MRW-SSA-Dvhop')
