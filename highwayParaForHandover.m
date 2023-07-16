%The function can be used to generate the corresponding parameters 
%of highway environments in the manuscript:

%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input: 
% L -> Number of APs 
% lengthOfWay -> Length of highway in simulation region 
% UE_lambda -> UE density
% numberOfLane -> Number of lane in highway

%Output:
% UEpositions -> User location coordinates, [K×1] vector
% APpositions -> AP location coordinates, [L×1] vector
% Garma -> Movement direction of UEs, [K×1] vector
  
function [UEpositions,APpositions,Garma] = highwayParaForHandover(L,lengthOfWay,UE_lambda,numberOfLane)
 
%The distance between AP with road in m 
AP_Road_Distance = 35; 
%The road width per line
road_Width_perLine = 4;  
%Inter-AP distance in m 
inter_AP_Dis = lengthOfWay/L;
%Minimum inter-distance of UE 
min_UE_Distance = 10;  
[UEpositions,Garma] = highwayUserPos(road_Width_perLine,lengthOfWay,UE_lambda,min_UE_Distance);   
%Calculate AP positions
road_Width = road_Width_perLine * numberOfLane;
% AP H- and V-axis
if L==1
    AP_Horizontal = [lengthOfWay/2];
else
    AP_Horizontal = [inter_AP_Dis/2 : inter_AP_Dis : L*inter_AP_Dis-inter_AP_Dis/2];
end
AP_Vertical = (AP_Road_Distance + road_Width) * ones(1,L);
APpositions = AP_Horizontal(:) + 1i*AP_Vertical(:);
% APpositions1 = APpositions-lengthOfWay;
% APpositions2 = APpositions+lengthOfWay;
% APpositions = [ APpositions1.' APpositions.' APpositions2.'];

 
