function [X,d]=MRW_Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData)
m=3;%将半径分为m级
%% 初始化节点间距离，跳数
clear Dall h
Dall = zeros(NodeAmount,NodeAmount);%距离矩阵
h = zeros(NodeAmount,NodeAmount);%跳数矩阵
for i = 1:NodeAmount
    for j = 1:NodeAmount
        Dall(i,j)=norm(data(2:3,i) - data(2:3,j),2); %计算距离
        if (Dall(i,j)<=R)&&(Dall(i,j)>0)
            h(i,j) = 1;
        elseif i == j
            h(i,j)=0;
        else
            h(i,j)=inf; %无效值
        end
    end
end
for i=1:BeaconAmount
        for j=1:NodeAmount
            for k=1:m
                if (Dall(i,j)<=k*R/4)&&(Dall(i,j)>(k-1)*R/4) 
                    h(i,j)=k/m;
                    break;
                end
            end
        end
end
%% 最短路径算法计算节点跳数
for k = 1:NodeAmount
    for i = 1:NodeAmount
        for j = 1:NodeAmount
            if h(i,k)+h(k,j)<h(i,j)
               h(i,j)=h(i,k)+h(k,j); 
            end
        end
    end
end
%% 计算每个信标节点的校正值(计算每个锚节点的平均每跳距离)
clear hBeacon DBeacon dhop
hBeacon = h(1:BeaconAmount,1:BeaconAmount);%跳数
DBeacon = Dall(1:BeaconAmount,1:BeaconAmount);%距离
for i = 1:BeaconAmount
   dhop(i) = sum(DBeacon(i,:))/sum(hBeacon(i,:)); 
end

%% 用跳数估计距离(计算每个未知节点到各锚节点的距离)
clear d hop1 hop a A B X1 X XWOA w hopsize
hop1 = h(1:BeaconAmount,(BeaconAmount+1):end);%未知节点到信标跳数，BeaconAmount行，UnAmount列 
%计算权值矩阵
for i=1:UnAmount
    w(:,i)=1./hop1(:,i)/sum(1./hop1(:,i)); %权值矩阵   
end
%根据权值计算未知节点的平均每跳距离
for i=1:UnAmount
      hopsize(1,i)=sum(w(:,i).*dhop(:,1));
end
for i=1:UnAmount
    hop=hopsize(1,i);  %hop为从最近信标获得的校正值
    d(:,i)=hop*hop1(:,i);    %BeaconAount行UnAmount列
end

%% 最小二乘法求未知节点坐标,构造A，B矩阵并求解
for i = 1:2
    for j=1:(BeaconAmount-1)
        a(i,j)=BeaconData(i,j)-BeaconData(i,BeaconAmount);
    end
end
A = -2*(a');
for j = 1:UnAmount
    for i = 1:(BeaconAmount-1)
        B(i,:)=d(i,j)^2-d(BeaconAmount,j)^2-BeaconData(1,i)^2+BeaconData(1,BeaconAmount)^2-BeaconData(2,i)^2+BeaconData(2,BeaconAmount)^2;
    end
     X1=inv(A'*A)*A'*B;
     X(:,j)=X1;
end
end