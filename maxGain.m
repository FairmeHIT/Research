%The function can be used to generate the corresponding parameters 
%of Algorithm1 in the manuscript:

%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
% K -> Number of UEs for the cell-free systems 
% L -> Number of APs  
% UE_flag -> %UE indicators that run to road boundaries, [K×1] vector 
% max_numAP -> Maximum number of serving APs, [K×L] matrix
% gainOverNoisedB -> Channel gains, [K×L] matrix 
% U_old -> The U matrix for the last time period, [K×L] matrix 
% set_AP_index_de_old -> AP subscript deleted in the last period, [K×L] matrix
  
%Output:
% U_new -> The U matrix of this time period, [K×L] matrix   
% set_AP_index_de_new -> AP index deleted in this period, [K×Ls] matrix  
% numHandover -> Number of handovers with LSFC-based method, [K×1] matrix  
% numPPHandover -> Number of PP handovers with LSFC-based method, [K×1] matrix  
   
function [U_new,set_AP_index_de_new,numHandover,numPPHandover] = ...
    maxGain(K,L,UE_flag,gainOverNoisedB,max_numAP, U_old, set_AP_index_de_old)

numPPHandover = zeros(K,1);
U_new = zeros(K,L);
numHandover = zeros(K,1); 
set_AP_index_de_new = zeros(K,max_numAP); 
for k = 1:K 
    if UE_flag(k,1) == 1
        U_new(k,:) = 0; 
        [~, I] = sort(gainOverNoisedB(k,:), 'descend'); 
        max_n_indices_temp = I(1:max_numAP); 
        U_new(k,max_n_indices_temp) = 1; 
    else 
        [~, I] = sort(gainOverNoisedB(k,:), 'descend'); 
        max_n_indices_new = I(1:max_numAP); 
        U_new(k,max_n_indices_new) = 1; 
        AP_index_new = max_n_indices_new; 
        [~,AP_index_old] = find(U_old(k,:));  
        %Added AP index
        setDiff_AP_index_add_new = setdiff(AP_index_new, AP_index_old); 
        numHandover(k) = length(setDiff_AP_index_add_new); 
        %Deleted AP index
        setDiff_AP_index_minus = setdiff(AP_index_old, AP_index_new); 
        set_AP_index_de_new(k,1:length(setDiff_AP_index_minus)) = setDiff_AP_index_minus; 
        %Whether the newly added items have been deleted last time
        Int_temp = intersect(setDiff_AP_index_add_new, set_AP_index_de_old(k,:)); 
        Int_temp(Int_temp==0)=[];
        numPPHandover(k) = length(Int_temp);  
    end 
end
end
     