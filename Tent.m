 %Tent shine upon
 function [x] = Tent(Max_iter)
 x(1)=rand; 
 a = 0.4;
 for i=1:Max_iter-1
     if x(i)<a
         x(i+1)=x(i)/a;
     end
     if x(i)>=a
         x(i+1)=(1-x(i)) / (1 - a);
     end
 end
 end