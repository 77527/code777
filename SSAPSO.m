
function [Best_score,Best_pos,curve]=SSAPSO(pop,Max_iter,lb,ub,dim,objfun)

%变量的范围
lb=lb*ones(1,dim); % 参数取值下界
ub=ub*ones(1,dim); % 参数取值上界

% 种群初始化
FS=initializationNew(pop,dim,ub,lb);  %%种群初始化
% 计算初始适应度值
fitness = zeros(1,pop);
for i = 1:pop %计算初始适应度
   fitness(i) =  objfun(FS(i,:));
end
 [~, index]= sort(fitness);%排序
GBestF = fitness(index(1));%全局最优适应度值
GBestX = FS(index(1),:);%全局最优位置
curve=zeros(1,Max_iter);
Fh =  FS(index(1),:); %核桃树上的松鼠
Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:);];%山核桃树上的松鼠
% 滑翔空气动力相关参数
row=1.204;V=5.25;S=0.0154;cd=0.6;CL=0.7;hg=1;sf=18;
Gc=1.9;Cr=0.5;	
D1=1/(2*row*V.^2*S*cd);L=1/(2*row*V.^2*S*CL);tanpi=D1/L;dg=hg/(tanpi*sf);
p1=0.15;p2=0.02; 
for t = 1: Max_iter
  disp(['当前迭代次数：',num2str(t)])
  pdp=(p1-p2)*(1-t/Max_iter)^5+p2;
  
    [~, index]= sort(fitness);%排序
    Fh = FS(index(1),:); %核桃树上的松鼠
    Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];%橡树上的松鼠
 for i = 1:3
       %% 橡树->山核桃树
      if (rand>= pdp)
       Fa(i,:) = Fa(i,:)+(dg*Gc*abs(Fh-Fa(i,:)));  
      else
           Fa(i,:) = Tent(dim).*(ub-lb)+lb;
       end
%        边界检查
      Temp = Fa(i,:);
      Temp(Temp>ub) = ub(Temp>ub);
      Temp(Temp<lb) = lb(Temp<lb);
      Fa(i,:) = Temp;
      FS(i+1,:)=Fa(i,:);
      fitness(i+1) = objfun(Fa(i,:));    
       [~, index]= sort(fitness);%排序
      Fh = FS(index(1),:);%山核桃树上的松鼠
      Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];%橡树上的松鼠
      if fitness(index(1))<GBestF
         GBestX = FS(index(1),:);
          GBestF = fitness(index(1));
      end
  end
   for i = 4:pop
      %% 普通树->橡树
    if(rand>pdp)
        FS(i,:)=FS(i,:)+(dg*Gc*abs(Fa(randi(3),:)-FS(i,:)));
      else
         FS(i,:)=Tent(dim).*(ub-lb)+lb;
     end
      %边界检查
      Temp = FS(i,:);
      Temp(Temp>ub) = ub(Temp>ub);
      Temp(Temp<lb) = lb(Temp<lb);
      FS(i,:) = Temp;
      fitness(i) = objfun(FS(i,:));
      [~, index]= sort(fitness);%排序
      Fh = FS(index(1),:);%山核桃树上的松鼠
      Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];%橡树树上的松鼠 
     if fitness(index(1))<GBestF
          GBestX= FS(index(1),:);
          GBestF = fitness(index(1));
     end 
      %% 普通树->山核桃树
      if(rand>pdp)
        FS(i,:)=FS(i,:)+(dg*Gc*abs(Fh-FS(i,:)));
      else
         FS(i,:)=Tent(dim).*(ub-lb)+lb;
       end
      %边界检查
      Temp = FS(i,:);
      Temp(Temp>ub) = ub(Temp>ub);
      Temp(Temp<lb) = lb(Temp<lb);
      FS(i,:) = Temp;
      fitness(i) = objfun(FS(i,:));
       [~, index]= sort(fitness);%排序
      Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];%橡树树上的松鼠
      Fh = FS(index(1),:);%山核桃树上的松鼠
     if fitness(index(1))<GBestF
          GBestX= FS(index(1),:);
          GBestF = fitness(index(1));
     end

   
      %% 季节变化条件
      Sc=sqrt(sum(sum(abs(Fh-Fa)).^2));
      Smin=(10*exp(-6))/(365).^(t/(Max_iter/2.5));
      if(Sc<Smin)
           FS(i,:)=lb+levy(1,dim,1.5).*(ub-lb);
           %边界检查
              Temp = FS(i,:);
              Temp(Temp>ub) = ub(Temp>ub);
              Temp(Temp<lb) = lb(Temp<lb);
              FS(i,:) = Temp;
              fitness(i) = objfun(FS(i,:));
              if fitness(i)<GBestF
                  GBestX= FS(i,:);
                  GBestF = fitness(i);
              end
      end   
   end
   %% 对FS进行高斯变异
     favg=mean(fitness);
   for i=1:pop 
       if fitness(i)<favg 
           k=normrnd(0,1);
           GSpop=FS(i,:).*(1+normrnd(0,1));   %%高斯变异（扰动）
           GS=objfun(GSpop);
           if GS<GBestF
               GBestX=GSpop;
               GBestF = GS;
           end
       end
   end
   
   %% PSO更新
   w = 0.9 - (0.9-0.4)*(t/Max_iter); % 惯性权重
   c1 = 2; % 个体学习因子
   c2 = 2; % 社会学习因子
   
   % 初始化速度(第一次迭代时)
   if t == 1
       V = zeros(pop, dim);
   end
   
   % PSO速度和位置更新
   for i = 1:pop
       % 速度更新
       V(i,:) = w*V(i,:) + c1*rand(1,dim).*(FS(i,:) - FS(i,:)) + ...
                c2*rand(1,dim).*(GBestX - FS(i,:));
       % 位置更新
       FS(i,:) = FS(i,:) + V(i,:);
       
       % 边界检查
       Temp = FS(i,:);
       Temp(Temp>ub) = ub(Temp>ub);
       Temp(Temp<lb) = lb(Temp<lb);
       FS(i,:) = Temp;
       
       % 更新适应度
       fitness(i) = objfun(FS(i,:));
       if fitness(i) < GBestF
           GBestX = FS(i,:);
           GBestF = fitness(i);
       end
   end
   
    curve(t) = GBestF;
end
Best_pos =GBestX;
Best_score = curve(end);
end