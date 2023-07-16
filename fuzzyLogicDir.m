%The function can be used to generate the corresponding parameters 
%of Algorithm1 in the manuscript:

%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
% K -> Number of UEs for the cell-free systems 
% U_old -> The U matrix for the last time period, [K×L] matrix 
% set_AP_index_de_old -> AP subscript deleted in the last period, [K×L] matrix
% max_numAP -> Maximum number of serving APs, [K×L] matrix
% channelGaindB -> Channel gains, [K×L] matrix
% Th_Gain -> HOM value, [K×1] vector
% Garma -> Movement direction of UEs, [K×1] vector
% Theta -> azimuth angle for (k,l) link, [K×L] matrix

 
%Output:
% U_new -> The U matrix of this time period, [K×L] matrix   
% set_AP_index_de_new -> AP index deleted in this period, [K×Ls] matrix  
% numHandover -> Number of handovers with Dir-based method, [K×1] matrix  
% numPPdirHandover -> Number of PP handovers with Dir-based method, [K×1] matrix  
 
function [U_new,set_AP_index_de_new,numHandover,numPPdirHandover] = ...
    fuzzyLogicDir(K,U_old,set_AP_index_de_old,max_numAP,channelGaindB,Th_HOM,Garma,Theta)
 
numHandover = zeros(K,1);
numPPdirHandover = zeros(K,1);
set_AP_index_de_new = zeros(K,max_numAP);
U_new = zeros(size(U_old));
for k = 1:K
    AP_index_old = find(U_old(k,:)==1);
    [old_min_gain, min_APindex_old_temp] = min(channelGaindB(k,AP_index_old));
    min_APindex_old = AP_index_old(min_APindex_old_temp); 
    U_temp_new = find(U_old(k,:)==0);
    [new_max_gain, max_APindex_new_temp] = max(channelGaindB(k,U_temp_new));
    max_APindex_new = U_temp_new(max_APindex_new_temp); 
    if Garma(k)==0 && rad2deg(Theta(k,max_APindex_new)) <= 90
         if  (new_max_gain - old_min_gain ) > Th_HOM(k)  
            numHandover(k) = 1;
            U_old(k,min_APindex_old) = 0;
            U_old(k,max_APindex_new) = 1;   
         end
    elseif Garma(k)~=0 && rad2deg(Theta(k,max_APindex_new)) >= 90
        if (new_max_gain - old_min_gain ) > Th_HOM(k)   
            numHandover(k) = 1;
            U_old(k,min_APindex_old) = 0;
            U_old(k,max_APindex_new) = 1;  
        end 
    end 
    U_new(k,:) = U_old(k,:);
    [~,AP_index_new] = find(U_new(k,:)==1);
    %Added AP index
    setDiff_AP_index_add_new = setdiff(AP_index_new, AP_index_old);  
    %Added AP number
    numHandover(k) = length(setDiff_AP_index_add_new); 
    %Deleted AP index
    setDiff_AP_index_minus = setdiff(AP_index_old, AP_index_new); 
    set_AP_index_de_new(k,1:length(setDiff_AP_index_minus)) = setDiff_AP_index_minus;
    %Whether the newly added items have been deleted last time
    Int_temp = intersect(setDiff_AP_index_add_new, set_AP_index_de_old(k,:)); 
    Int_temp(Int_temp==0)=[];
    numPPdirHandover(k) = length(Int_temp); 
end
end


