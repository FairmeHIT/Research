%The function can be used to calculate the analytical result of SE_k[n] 
%according to Theorem 1 in the manuscript:
 
%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
% K -> Number of UEs for the cell-free systems
% M -> Antenna number of AP
% L -> Number of APs    
% U -> Serving AP indicator matrix, U(k,l)=1 means that l is the serving AP of k, generated by (54), [K×L] matrix 
% P -> Transmit power of UE in W, [K×1] vector
% Hhat_LMMSE -> Estimate channel matrix generated by <channelEstimateLMMSE.m>, [M*L×N×K] matrix 
% H -> Channel matrix generated by Eq.(9), [M*L×N×K] matrix  
% Omega -> Spatial correlation matrix generated by <scfGenerate.m>, [M×M×L×K] matrix 
% R -> Spatial correlation matrix generated by <scfGenerate.m>, [M×M×L×K] matrix 
% B_LMMSE -> Estimate covariance matrix with LMMSE estimator generated by <channelEstimateLMMSE.m>, [M*L×N×K] matrix
% ACF_2_d -> ACF_2 for S=[T+1:C] generated by <acfGenerate.m>, [K×L×T] vector
% ACF_2_t -> ACF_2 for S=[1+T] generated by <acfGenerate.m>, [K×L×T] vector
% a_h -> Wave vector generated by Eq.(2), [M×L×K] vector
% J -> Intermediate variable given by <channelEstimateLMMSE>, [M×M×L×K] matrix 
% n -> n-th symbol of transmission phase

%Output:
% SE_analytical_kn_Theorem2, LMMSE estimator, n-th symbol, [K×1] vector   

function [SE_analytical_kn_Theorem2] = analyticalSEwithLMMSE(K,M,L,U,P,Hhat_LMMSE,Omega,R,B_LMMSE,ACF_2_d,ACF_2_t,a_h,J,n)

%Prepare to store
Hhat_Aging = zeros(size(Hhat_LMMSE));  
B_Aging = zeros(size(B_LMMSE));  
for k = 1:K 
    for l = 1:L 
        %intermediate variable
        Hhat_Aging((l-1)*M_AP+1:l*M_AP,:,k) = ACF_2_d(k,l,n)*Hhat_LMMSE((l-1)*M+1:l*M,:,k);
        B_Aging(:,:,l,k) = (abs(ACF_2_d(k,l,n)))^2 * B_LMMSE(:,:,l,k); 
    end
end
 
%Prepare to store the terms that appear in SEs
los_u_2 = zeros(K,L,K);
los_u_1 = zeros(K,L,K);
Nois = zeros(L,K);   
%Prepare to store the SE
SE_analytical_kn_Theorem2 = zeros(K,1);  
  
%Compute the closed-form expression according to Theorem 2
%Go through each AP
for l = 1:L
    %Extract which UEs are served by the AP l
    UEs_index = find(U(:,l)==1); 
    %Go through all UEs served by the AP l
    for UE_i = 1:length(UEs_index) 
        %Extract UE index
        k = UEs_index(UE_i); 
        %Noise scaling  
        Nois(l,k) =  (trace(B_Aging(:,:,l,k))); 
        %Compute the received pilot signal covariance matrix Psi
        Psi = (P*sum(Omega(:,:,l,pilotIndex(k)==pilotIndex),4) + eyeN); 
        %Go through all UEs
        for i = 1:K 
            %Eq.(70)
            los_u_1(i,l,k) =  (trace(B_Aging(:,:,l,k)*Omega(:,:,l,i)));  
            %If UE i shares the same pilot with UE k
            if pilotIndex(k) == pilotIndex(i) 
                %Eq.(43) but using ACF_2 
                u_l = P * ACF_2_d(i,l,n)*ACF_2_t(k,l,pilotIndex(k))*ACF_2_t(i,l,pilotIndex(k))'; 
                %Eq.(86)
                los_u_2(i,l,k) =   (u_l) *  (trace(Omega(:,:,l,k)*J(:,:,l,i))); 
                %a1 and a2 for Eq.(75)
                %Eq.(78)
                a1 = trace(J(:,:,l,k)*R(:,:,l,i)*J(:,:,l,k)'*(Psi-P*Omega(:,:,l,i))) + ...
                    a_h(:,l,i)'*J(:,:,l,k)'*(Psi-P*Omega(:,:,l,i))*J(:,:,l,k)*a_h(:,l,i);
                %Eq.(80)
                a2 = P*abs(ACF_2_t(i,l,pilotIndex(i)))^2*( abs(trace(R(:,:,l,i)*J(:,:,l,k)))^2 + ... 
                    2*real(trace(R(:,:,l,i)*J(:,:,l,k))*a_h(:,l,i)'*J(:,:,l,k)'*a_h(:,l,i))) + ...
                    P*trace(Omega(:,:,l,i)*J(:,:,l,k)*Omega(:,:,l,i)*J(:,:,l,k)') ;  
                %Eq.(85)
                r_1 = abs(ACF_2_t(k,l,pilotIndex(k)))^2*P*(a1 + a2) ; 
                %for Eq.(72)
                los_u_1(i,l,k) = los_u_1(i,l,k) + abs(ACF_2_d(i,l,n))^2 *( (r_1) - ...
                    (trace(B_Aging(:,:,l,k)*Omega(:,:,l,i))))  ;
            end 
        end 
    end 
end 
% Compute the SEs  
for k = 1:K 
    %Determine the set of serving APs for UE k
    seAPs = find(U(k,:)==1); 
    %for C_ki^{los-u}
    vec_temp = vec(sqrt(P)*los_u_2(k,seAPs,k)); 
    Mat_temp =  los_u_2(:,seAPs,k).'*conj(los_u_2(:,seAPs,k));
    denom_matrix = P*(diag(sum(los_u_1(:,seAPs,k),1))+Mat_temp-diag(diag(Mat_temp)))...
        -vec_temp*vec_temp'+diag(Nois(seAPs,k)); 
    %%Compute the LSFD vector 
    a_k = denom_matrix\vec_temp; 
    %Compute the SE for symbol n 
    SE_analytical_kn_Theorem2(k) = real(log2(1+abs(a_k'*vec_temp)^2/(a_k'*denom_matrix*a_k))); 
end

 
 