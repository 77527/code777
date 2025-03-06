function [X,d]=Dvhop(BeaconAmount,UnAmount,NodeAmount,R,data,BeaconData)
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
clear hBeacon DBeacon dhop DUN dhopUN
hBeacon = h(1:BeaconAmount,1:BeaconAmount);%����
DBeacon = Dall(1:BeaconAmount,1:BeaconAmount);%����
for i = 1:BeaconAmount
   dhop(i) = sum(DBeacon(i,:))/sum(hBeacon(i,:)); 
end
%% ����ÿ��δ֪�ڵ��У��ֵ
DUN =  Dall(1:BeaconAmount,(BeaconAmount+1):NodeAmount);    %BeaconAmount��UNAmount��  �ű�ڵ㵽δ֪�ڵ�ľ���
for i=1:BeaconAmount
    for j = 1:UnAmount
        if min(DUN(:,j))==DUN(i,j)
            dhopUN(j) = dhop(i);%δ֪�ڵ��������ű���У��ֵ��
        end                     %�����δ֪�ڵ������ê�ڵ��ƽ��ÿ��������Ϊ�ýڵ��ÿ������        
    end
end

%% ���������ƾ���(����ÿ��δ֪�ڵ㵽��ê�ڵ�ľ���)
clear d hop1 hop a A B X1 X
hop1 = h(1:BeaconAmount,(BeaconAmount+1):end);%δ֪�ڵ㵽�ű�������BeaconAmount�У�UnAmount�� 
for i=1:UnAmount
    hop=dhopUN(i);  %hopΪ������ű��õ�У��ֵ
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