clear;
close all;
clc;

NodeAmountArr=[100 150 200 250 300 350 400] ; % �ܽڵ���������
BorderLength = 100; %����߽緶Χ
times = 1; %���д���

for run=1:times
for index = 1:length(NodeAmountArr)
NodeAmount=NodeAmountArr(index);
%������Χ��������ɽڵ㣬���ܽڵ���NodeAmount������
AreaC = BorderLength.*rand(2,NodeAmount);%[x1,...,xn;y1,...,yn;];
BeaconAmount = NodeAmount * 0.2; %�ű�ڵ�����ê�ڵ㣩
UnAmount = NodeAmount - BeaconAmount; %δ֪�ڵ���
R = 30; %ͨ�ž���
%Ϊÿ���������ţ����1��2��3�����ڵ�1��
data = [(1:NodeAmount);AreaC];
%�ű�������Ϣ
BeaconData = data(2:3,1:BeaconAmount);%��ȡ2��3�д�ŵ�����
UnKnownData = data(2:3,BeaconAmount+1:end);%��ȡʣ�µ�����Ϊδ֪�ڵ�����

%% �㷨1��ԭʼDvhop
[X,d]=Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData);

%% �㷨2����ͨ�Ű뾶�������Ȩ��Dvhop
[X_w,d_w]=MRW_Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData);

%% �㷨3��SSAPSO-Dvhop
%% ���������㷨����
XSSA = [];
for i = 1:UnAmount
    dim = 2; %dim�Ǳ����������������ά�ȣ�
    pop = 40;%��Ⱥ����
    MaxIter = 50;%����������
    ub = BorderLength; %�������Χ
    lb = 0;%������С��Χ
    fobj = @(x) fun(x,UnAmount,BeaconAmount,d,BeaconData,i);%��Ӧ�Ⱥ������� 
    [Best_score,Best_pos,SSA_curve(i,:)]=SSAPSO(pop,MaxIter,lb,ub,dim,fobj); %��ʼ�Ż�
    XSSA = [XSSA;Best_pos];
end
XSSA =XSSA';

%% �㷨4��MRW-SSAPSO-Dvhop
%% ���������Ż��㷨����
XWSSA = [];
for i = 1:UnAmount
    dim = 2; %dim�Ǳ����������������ά�ȣ�
    pop = 20;%��Ⱥ����
    MaxIter = 20;%����������
    ub = BorderLength; %�������Χ
    lb = 0;%������С��Χ
    fobj = @(x) fun(x,UnAmount,BeaconAmount,d_w,BeaconData,i);%��Ӧ�Ⱥ������ã��õ���d_w
    [Best_score_WSSA,Best_pos_WSSA,WSSA_curve(i,:)]=SSAPSO(pop,MaxIter,lb,ub,dim,fobj); %��ʼ�Ż�
    XWSSA = [XWSSA;Best_pos_WSSA];
end
XWSSA =XWSSA';

%�������
clear error errorSSA error_WSSA
for i=1:UnAmount
    error(1,i)=((X(1,i)-UnKnownData(1,i))^2+(X(2,i)-UnKnownData(2,i))^2)^0.5;
    errorSSA(1,i) = ((XSSA(1,i)-UnKnownData(1,i))^2+(XSSA(2,i)-UnKnownData(2,i))^2)^0.5;
    error_WSSA(1,i) = ((XWSSA(1,i)-UnKnownData(1,i))^2+(XWSSA(2,i)-UnKnownData(2,i))^2)^0.5;
end

%����Dvhop��������
Accuracy(run,index)=sum(error)/(UnAmount*R);%��һ����λ���
%SSA_Dvhop��������
AccuracySSA(run,index)=sum(errorSSA)/(UnAmount*R);
%MRW-ISSA-hop��������
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
% legend('ԭʼDvhop','SSA-Dvhop','MRW-SSA-Dvhop')