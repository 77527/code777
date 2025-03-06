%% The fitness function
function fitness =  fun(x,UnAmount,BeaconAmount,d,BeaconData,index)
% x = reshape(x,[2,UnAmount]);
TempValue = 0;

for i = 1:(BeaconAmount-1)
       TempValue =TempValue +  abs(sqrt((x(1)-BeaconData(1,i))^2 + (x(2)-BeaconData(2,i))^2) - d(i,index));
end

fitness = TempValue;
end