%The function can be used to generate the corresponding parameters 
%of highway environments in the manuscript:

%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

% Input: 
% width_perLane -> Width of each lane  
% lengthOfWay -> Length of highway in simulation region 
% UE_lambda -> UE density in /m/lane
% min_Distance -> Minimum spaceing between UE

% Output:
% UEpositions -> User location coordinates, [K×1] vector 
% Garma -> Movement direction of UEs, [K×1] vector

function [UEpositions,Garma] = highwayUserPos(width_perLane,lengthOfWay,UE_lambda,min_Distance)
 
% For Line 01
massUEpoint = lengthOfWay * UE_lambda;
numLinePoint = poissrnd(massUEpoint);
line1_x=zeros(numLinePoint,1);
i = 1;
while i<=numLinePoint
    line1_x(i) = rand(1,1)*lengthOfWay;
    if(i>=2)
        for j=1:i-1
            if(abs(line1_x(i)-line1_x(j)) < min_Distance)
                i = i - 1;
                break;
            end
        end
    end
    i = i + 1;
end
line1_y = (width_perLane/2 + (1-1)*width_perLane)*ones(numLinePoint,1);
Garma1 = 0*ones(length(line1_y),1); 
% For Line 02
massUEpoint = lengthOfWay * UE_lambda;
numLinePoint = poissrnd(massUEpoint);
line2_x=zeros(numLinePoint,1);
i = 1;
while i<=numLinePoint
    line2_x(i) = rand(1,1)*lengthOfWay;
    if(i>=2)
        for j=1:i-1
            if(abs(line2_x(i)-line2_x(j)) < min_Distance)
                i = i - 1;
                break;
            end
        end
    end
    i = i + 1;
end
line2_y = (width_perLane/2 + (2-1)*width_perLane)*ones(numLinePoint,1);
Garma2 = 0*ones(length(line2_y),1);
% For Line 03
massUEpoint = lengthOfWay * UE_lambda;
numLinePoint = poissrnd(massUEpoint);
line3_x=zeros(numLinePoint,1);
i = 1;
while i<=numLinePoint
    line3_x(i) = rand(1,1)*lengthOfWay;
    if(i>=2)
        for j=1:i-1
            if(abs(line3_x(i)-line3_x(j)) < min_Distance)
                i = i - 1;
                break;
            end
        end
    end
    i = i + 1;
end
line3_y = (width_perLane/2 + (3-1)*width_perLane)*ones(numLinePoint,1);
Garma3 = 0*ones(length(line3_y),1);
% For Line 04
massUEpoint = lengthOfWay * UE_lambda;
numLinePoint = poissrnd(massUEpoint);
line4_x=zeros(numLinePoint,1);
i = 1;
while i<=numLinePoint
    line4_x(i) = rand(1,1)*lengthOfWay;
    if(i>=2)
        for j=1:i-1
            if(abs(line4_x(i)-line4_x(j)) < min_Distance)
                i = i - 1;
                break;
            end
        end
    end
    i = i + 1;
end
line4_y = (width_perLane/2 + (4-1)*width_perLane)*ones(numLinePoint,1);
Garma4 = pi*ones(length(line4_y),1);
% For Line 05
massUEpoint = lengthOfWay * UE_lambda;
numLinePoint = poissrnd(massUEpoint);
line5_x=zeros(numLinePoint,1);
i = 1;
while i<=numLinePoint
    line5_x(i) = rand(1,1)*lengthOfWay;
    if(i>=2)
        for j=1:i-1
            if(abs(line5_x(i)-line5_x(j)) < min_Distance)
                i = i - 1;
                break;
            end
        end
    end
    i = i + 1;
end
line5_y = (width_perLane/2 + (5-1)*width_perLane)*ones(numLinePoint,1);
Garma5 = pi*ones(length(line5_y),1);
% For Line 06
massUEpoint = lengthOfWay * UE_lambda;
numLinePoint = poissrnd(massUEpoint);
line6_x=zeros(numLinePoint,1);
i = 1;
while i<=numLinePoint
    line6_x(i) = rand(1,1)*lengthOfWay;
    if(i>=2)
        for j=1:i-1
            if(abs(line6_x(i)-line6_x(j)) < min_Distance)
                i = i - 1;
                break;
            end
        end
    end
    i = i + 1;
end
line6_y = (width_perLane/2 + (6-1)*width_perLane)*ones(numLinePoint,1);
Garma6 = pi*ones(length(line6_y),1);
Line_x = [line1_x;line2_x;line3_x;line4_x;line5_x;line6_x];
Line_y = [line1_y;line2_y;line3_y;line4_y;line5_y;line6_y];
Garma = [Garma1;Garma2;Garma3;Garma4;Garma5;Garma6];  
UEpositions = Line_x + 1i*Line_y;
end

 