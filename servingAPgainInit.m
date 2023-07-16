%The function can be used to design the initial U matrix in the manuscript:
 
%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.
 
%Input:
% K -> Number of UEs for the cell-free systems 
% L -> Number of APs 
% gainOverNoisedB -> Channel gains, [K×L] matrix 
% max_numAP -> Maximum number of serving APs, [K×L] matrix
 
%Output:
% U_new -> The U matrix of this time period, [K×L] matrix   
% set_AP_index_de_new -> AP index added in this period, [K×Ls] matrix   
 
function [U_AP_kl_new,set_AP_index_de] = servingAPgainInit(K,L,gainOverNoisedB,max_numAP)
%Prepare to store results
U_AP_kl_new = zeros(K,L); 
max_n_indices_old = [];   
set_AP_index_de = zeros(K,max_numAP); 
for k = 1:K
    [~, I] = sort(gainOverNoisedB(k,:), 'descend');
    max_n_indices_new = I(1:max_numAP);
    U_AP_kl_new(k,max_n_indices_new) = 1;
    index_temp = setdiff(max_n_indices_old, max_n_indices_new);
    %Added AP index
    if isempty(index_temp)
        set_AP_index_de(k,:) = 0; %Added AP index
    else
        set_AP_index_de(k,:) = index_temp;
    end 
end
end