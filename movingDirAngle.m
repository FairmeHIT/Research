%The function can be used to generate the corresponding parameters 
%of Algorithm1 in the manuscript:

%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
% ueCoordinate1 -> P_k(x_i,y_i) 
% ueCoordinate1 -> P_k(x_i-1,y_i-1) 
% apCoordinate -> O_l(x_0,y_0)
 
%Output:
% vartheta -> vartheta value in [deg] 

function vartheta = movingDirAngle(ueCoordinate1,ueCoordinate2,apCoordinate)

if ueCoordinate1(1)<ueCoordinate2(1)
        beta=atand((ueCoordinate2(2)-ueCoordinate1(2))/(ueCoordinate2(1)-ueCoordinate1(1)));
else
        beta=180+atand((ueCoordinate2(2)-ueCoordinate1(2))/(ueCoordinate2(1)-ueCoordinate1(1)));
end
        d=abs(tan(beta)*(ueCoordinate2(1)-apCoordinate(1))+(ueCoordinate2(2)-apCoordinate(2)))./(sqrt((tan(beta)).^2+1));
if ueCoordinate2(1)<=apCoordinate(1)
        vartheta=asind(d./(sqrt((ueCoordinate2(1)-apCoordinate(1)).^2+(ueCoordinate2(2)-apCoordinate(2)).^2)));
else 
        vartheta=180-asind(d./(sqrt((ueCoordinate2(1)-apCoordinate(1)).^2+(ueCoordinate2(2)-apCoordinate(2)).^2)));
end 
end
 