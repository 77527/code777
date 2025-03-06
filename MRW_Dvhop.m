function [X,d]=MRW_Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData)
m=3;%���뾶��Ϊm��
%% ��ʼ���ڵ����룬����
clear Dall h
Dall = zeros(NodeAmount,NodeAmount);%�������
h = zeros(NodeAmount,NodeAmount);%��������
for i = 1:NodeAmount
    for j = 1:NodeAmount
        Dall(i,j)=norm(data(2:3,i) - data(2:3,j),2); %�������
        if (Dall(i,j)<=R)&&(Dall(i,j)>0)
            h(i,j) = 1;
        elseif i == j
            h(i,j)=0;
        else
            h(i,j)=inf; %��Чֵ
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
%% ���·���㷨����ڵ�����
for k = 1:NodeAmount
    for i = 1:NodeAmount
        for j = 1:NodeAmount
            if h(i,k)+h(k,j)<h(i,j)
               h(i,j)=h(i,k)+h(k,j); 
            end
        end
    end
end
%% ����ÿ���ű�ڵ��У��ֵ(����ÿ��ê�ڵ��ƽ��ÿ������)
clear hBeacon DBeacon dhop
hBeacon = h(1:BeaconAmount,1:BeaconAmount);%����
DBeacon = Dall(1:BeaconAmount,1:BeaconAmount);%����
for i = 1:BeaconAmount
   dhop(i) = sum(DBeacon(i,:))/sum(hBeacon(i,:)); 
end

%% ���������ƾ���(����ÿ��δ֪�ڵ㵽��ê�ڵ�ľ���)
clear d hop1 hop a A B X1 X XWOA w hopsize
hop1 = h(1:BeaconAmount,(BeaconAmount+1):end);%δ֪�ڵ㵽�ű�������BeaconAmount�У�UnAmount�� 
%����Ȩֵ����
for i=1:UnAmount
    w(:,i)=1./hop1(:,i)/sum(1./hop1(:,i)); %Ȩֵ����   
end
%����Ȩֵ����δ֪�ڵ��ƽ��ÿ������
for i=1:UnAmount
      hopsize(1,i)=sum(w(:,i).*dhop(:,1));
end
for i=1:UnAmount
    hop=hopsize(1,i);  %hopΪ������ű��õ�У��ֵ
    d(:,i)=hop*hop1(:,i);    %BeaconAount��UnAmount��
end

%% ��С���˷���δ֪�ڵ�����,����A��B�������
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