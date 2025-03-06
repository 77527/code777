
function [Best_score,Best_pos,curve]=SSAPSO(pop,Max_iter,lb,ub,dim,objfun)

%range of variables 
lb=lb*ones(1,dim); 
ub=ub*ones(1,dim); 

% Tent initializes the population
FS=initializationNew(pop,dim,ub,lb);  
% The initial fitness value is calculated
fitness = zeros(1,pop);
for i = 1:pop 
   fitness(i) =  objfun(FS(i,:));
end
 [~, index]= sort(fitness);
GBestF = fitness(index(1));
GBestX = FS(index(1),:);
curve=zeros(1,Max_iter);
Fh =  FS(index(1),:); 
Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:);];
% Gliding aerodynamic-related parameters
row=1.204;V=5.25;S=0.0154;cd=0.6;CL=0.7;hg=1;sf=18;
Gc=1.9;Cr=0.5;	
D1=1/(2*row*V.^2*S*cd);L=1/(2*row*V.^2*S*CL);tanpi=D1/L;dg=hg/(tanpi*sf);
p1=0.15;p2=0.02; 
for t = 1: Max_iter
  disp(['Number of current iterations：',num2str(t)])
  pdp=(p1-p2)*(1-t/Max_iter)^5+p2;
  
    [~, index]= sort(fitness);
    Fh = FS(index(1),:); 
    Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];
 for i = 1:3
       %%oak->hickory 
      if (rand>= pdp)
       Fa(i,:) = Fa(i,:)+(dg*Gc*abs(Fh-Fa(i,:)));  
      else
           Fa(i,:) = Tent(dim).*(ub-lb)+lb;
       end
%      bounds checking 
      Temp = Fa(i,:);
      Temp(Temp>ub) = ub(Temp>ub);
      Temp(Temp<lb) = lb(Temp<lb);
      Fa(i,:) = Temp;
      FS(i+1,:)=Fa(i,:);
      fitness(i+1) = objfun(Fa(i,:));    
       [~, index]= sort(fitness);
      Fh = FS(index(1),:);
      Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];
      if fitness(index(1))<GBestF
         GBestX = FS(index(1),:);
          GBestF = fitness(index(1));
      end
  end
   for i = 4:pop
      %% Ordinary tree->oak
    if(rand>pdp)
        FS(i,:)=FS(i,:)+(dg*Gc*abs(Fa(randi(3),:)-FS(i,:)));
      else
         FS(i,:)=Tent(dim).*(ub-lb)+lb;
     end
      %bounds checking
      Temp = FS(i,:);
      Temp(Temp>ub) = ub(Temp>ub);
      Temp(Temp<lb) = lb(Temp<lb);
      FS(i,:) = Temp;
      fitness(i) = objfun(FS(i,:));
      [~, index]= sort(fitness);
      Fh = FS(index(1),:);
      Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];
     if fitness(index(1))<GBestF
          GBestX= FS(index(1),:);
          GBestF = fitness(index(1));
     end 
      %%Ordinary tree->hickory 
      if(rand>pdp)
        FS(i,:)=FS(i,:)+(dg*Gc*abs(Fh-FS(i,:)));
      else
         FS(i,:)=Tent(dim).*(ub-lb)+lb;
       end
      %bounds checking
      Temp = FS(i,:);
      Temp(Temp>ub) = ub(Temp>ub);
      Temp(Temp<lb) = lb(Temp<lb);
      FS(i,:) = Temp;
      fitness(i) = objfun(FS(i,:));
       [~, index]= sort(fitness);
      Fa = [FS(index(2),:);FS(index(3),:);FS(index(4),:)];
      Fh = FS(index(1),:);
     if fitness(index(1))<GBestF
          GBestX= FS(index(1),:);
          GBestF = fitness(index(1));
     end

   
      %% Section change conditions
      Sc=sqrt(sum(sum(abs(Fh-Fa)).^2));
      Smin=(10*exp(-6))/(365).^(t/(Max_iter/2.5));
      if(Sc<Smin)
           FS(i,:)=lb+levy(1,dim,1.5).*(ub-lb);
           %bounds checking
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
   %% FS performed the Gaussian variation
     favg=mean(fitness);
   for i=1:pop 
       if fitness(i)<favg 
           k=normrnd(0,1);
           GSpop=FS(i,:).*(1+normrnd(0,1));  
           GS=objfun(GSpop);
           if GS<GBestF
               GBestX=GSpop;
               GBestF = GS;
           end
       end
   end
   
   %% PSO update
   w = 0.9 - (0.9-0.4)*(t/Max_iter); 
   c1 = 2; 
   c2 = 2; 
   
   % Initialization speed
   if t == 1
       V = zeros(pop, dim);
   end
   
   % PSO speed and location updates
   for i = 1:pop
       V(i,:) = w*V(i,:) + c1*rand(1,dim).*(FS(i,:) - FS(i,:)) + ...
                c2*rand(1,dim).*(GBestX - FS(i,:));
       FS(i,:) = FS(i,:) + V(i,:);
       
       % bounds checking
       Temp = FS(i,:);
       Temp(Temp>ub) = ub(Temp>ub);
       Temp(Temp<lb) = lb(Temp<lb);
       FS(i,:) = Temp;
       
       % Update fitness
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