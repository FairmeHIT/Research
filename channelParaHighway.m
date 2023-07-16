%The function can be used to generate the propagation parameters in the manuscript:

%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input: 
% L -> Number of APs 
% UEpositions -> User location coordinates, [K×1] vector
% APpositions -> AP location coordinates, [L×1] vector
 
%Output:
% probLOS_allLoS -> LoS probability for highway, [K×L] matrix 
% ricianFactor_allLoS -> Rician factor for highway, [K×L] matrix  
% channelGaindB_allLoS -> Channel gain in dB, [K×L] matrix  
% probLOS_allLoS -> LoS probability for Rayleigh fading, [K×L] matrix 
% ricianFactor_allLoS -> Rician factor for Rayleigh fading, [K×L] matrix  
% channelGaindB_allLoS -> Channel gain Rayleigh fading, [K×L] matrix  
 
function [probLOS_allLoS,ricianFactor_allLoS,channelGaindB_allLoS,...
    probLOS_nonLoS,ricianFactor_nonLoS,channelGaindB_nonLoS,Theta] = channelParaHighway(L,UEpositions,APpositions)

%Communication bandwidth
B = 20e6; 
%Noise figure (in dB)
noiseFigure = 7;
%AP height (in m)
AP_h = 5; 
%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure; 
%Calculate the UE number
K = size(UEpositions,1);  
gainGaindB_kl_allLoS = zeros(K,L); 
gainGaindB_kl_nonLoS = zeros(K,L);   
Theta = zeros(K,L);  
%Standard deviation of shadow fading in dB
sigma_sf_NLOS=8; %for NLOS
sigma_sf_LOS=3;   %for LOS  
probLOS_allLoS=ones(K,L); 
probLOS_nonLoS=zeros(K,L);
ricianFactor_nonLoS=zeros(K,L);
ricianFactor = zeros(K,L);  
randn_temp = randn(K,L);
 
for k = 1:K
    %2D distance between UE and AP
    [distanceto_k] =  abs(APpositions - UEpositions(k));  
    %3D distance
    distances_k_3D = sqrt(AP_h^2 + distanceto_k.^2); 
    ricianFactor(k,:)=db2pow(13-0.03*distances_k_3D);   
    for l = 1:L 
        % for allLoS in highway environment
        gainGaindB_kl_allLoS(k,l) = -30.18 - 26*log10(distances_k_3D(l)); 
        shadowing = sigma_sf_LOS*randn_temp(k,l);
        channelGainShadowing = gainGaindB_kl_allLoS(k,l) + shadowing; 
        gainGaindB_kl_allLoS(k,l) = channelGainShadowing;    
        % for nonLoS in highway environment
        gainGaindB_kl_nonLoS(k,l) =  -34.53 - 38*log10(distances_k_3D(l)); 
        shadowing = sigma_sf_NLOS*randn_temp(k,l);
        channelGainShadowing = gainGaindB_kl_nonLoS(k,l) + shadowing; 
        gainGaindB_kl_nonLoS(k,l) = channelGainShadowing;    
        %Calculate Theta
        Theta(k,l) = angle(APpositions(l) - UEpositions(k));   
    end 
end 

ricianFactor_allLoS = ricianFactor;     
channelGaindB_allLoS = gainGaindB_kl_allLoS-noiseVariancedBm;           
channelGaindB_nonLoS = gainGaindB_kl_nonLoS-noiseVariancedBm;   
end
 