%The script can be used to count the number of handovers 
%according to the soft handover algorithm in the manuscript:
 
%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
%None
  
%Output:
% Num_handover_kN_fuzzyBased_Dir -> Number of handovers with fuzzy&Dir-based method for different UE densities
% Num_PP_handover_kN_fuzzyBased_Dir -> Number of PP handovers with fuzzy&Dir-based method for different UE densities
% Num_handover_kN_fuzzyBased_nonDir -> Number of handovers with fuzzy-based method for different UE densities
% Num_PP_handover_kN_fuzzyBased_nonDir -> Number of PP handovers with fuzzy-based method for different UE densities
% Num_handover_kN_LSFC -> Number of handovers with LSFC-based method for different UE densities
% Num_PP_handover_kN_LSFC -> Number of PP handovers with LSFC-based method for different UE densities

clc,clear all,close all 
warning('off');
C = 310;
%symbol length
T_s=10^(-5); 
%pilot number
T = 10;
DataSym = C - T; 
T_t = 50; %in s  
%Upgrade frequency
update_Freq = 1; %in s
%n' range
n_p_r = [0:update_Freq:T_t]; 
%parameters with varying 
Rician_range = [1];   
max_num_APserving_range = [5]; %for all AP in this network 
%freeway length (simulation region)
lengthOfWay = 3464;  
%Number of lane
road_Num = 6;
%namuda in [m] wavelenght of 2GHz
namuda = 0.15;  
%UE density range, in UE/lane/m
UE_lambda_r = [0.001:0.001:0.01];
%Number of AP 
L = 20;  
N = 100;

% Go through all realizations
for n = 1:N
    % Go through all UE density 
    for UE_lambda_i = 1:length(UE_lambda_r) 
        UE_lambda = UE_lambda_r(UE_lambda_i); 
        [UEpositions,APpositions,Garma] = highwayParaForHandover(L,lengthOfWay,UE_lambda,road_Num);   
        %Number of UE
        K = length(UEpositions);  
        %The moving speed of UE is set randomly and uniformly distributed in the range of [70~140] km/h for freeway scenarios  
        v_range = [70 140]/3.6;%in m/s 
        v_r = rand(K,1)*(max(v_range)-min(v_range)) + min(v_range);  
        %The spread of AoA is chosen randomly and uniformly distributed in the range of [10^{\circ}~100^{\circ}]
        UE_angleSpread = [100 10];  
        ka_range_UE = (1./deg2rad(UE_angleSpread)).^2;  
        kappa_T_r = rand(K,1)*(max(ka_range_UE)-min(ka_range_UE)) + min(ka_range_UE);  

        for serving_num_AP_i = 1:length(max_num_APserving_range)  
            max_num_APserving = max_num_APserving_range(serving_num_AP_i); 
            for n_i = 1:length(n_p_r)   
                  disp(['  [n_i = ' num2str(n_i) ', All = ' num2str(length(n_p_r)) ']'... 
                        '  [serving_num_AP_i = ' num2str(serving_num_AP_i) ', All = ' num2str(length(max_num_APserving_range)) ']'...
                        '  [UE_lambda_i = ' num2str(UE_lambda_i) ', All = ' num2str(length(UE_lambda_r)) ']'...  
                        '  [n = ' num2str(n) ', All = ' num2str(N) ']']); 
                %Initialization indicator
                n_yipi = n_p_r(n_i);
                %Prepare to store the UE that moving to the road boundary
                UE_flag = zeros(K,1);  
                %Record the actual location of the UEs
                for k=1:K  
                    v = v_r(k);  
                    if Garma(k) == 0 
                        UEpositions(k,1) = (real(UEpositions(k,1)) + v*update_Freq) + 1i*imag(UEpositions(k,1)); 
                    else  
                        UEpositions(k,1) = (real(UEpositions(k,1)) - v*update_Freq) + 1i*imag(UEpositions(k,1));
                    end  
                    if real(UEpositions(k,1)) > lengthOfWay  
                        UE_flag(k,1) = 1; %Record the user to the boundary
                        UEpositions(k,1) = 0 + 1i*imag(UEpositions(k,1));
                    elseif real(UEpositions(k,1)) < 0  
                        UE_flag(k,1) = 1; %Record the user to the boundary
                        UEpositions(k,1) = lengthOfWay + 1i*imag(UEpositions(k,1));
                    end 
                end  
                [probLOS_allLoS,ricianFactor_allLoS,channelGaindB_allLoS,...
                probLOS_nonLoS,ricianFactor_nonLoS,channelGaindB_nonLoS,Theta] = channelParaHighway(L,UEpositions,APpositions); 
                % Go through all fading (Rician or Rayleigh)
                for rician_i = 1:length(Rician_range) 
                    rician_flag = Rician_range(rician_i); 
                    if rician_flag == 1
                        fis=readfis('handoverFuzzyRician'); 
                        probLOS = probLOS_allLoS;
                        ricianFactor = ricianFactor_allLoS;
                        channelGaindB = channelGaindB_allLoS;  
                    else
                        fis=readfis('handoverFuzzyRayleigh'); 
                        probLOS = probLOS_nonLoS;
                        ricianFactor = ricianFactor_nonLoS;
                        channelGaindB = channelGaindB_nonLoS;  
                    end 
                    TCC_Matrix_Data_d =  zeros(K,L,DataSym); 
                    TCC_Matrix_Pilot_t = zeros(K,L,T);   
                    % Go through all UE for ACF calculation 
                    for k = 1:K  
                        ka_UE = kappa_T_r(k); 
                        fD = v_r(k)/namuda;
                        fDTs = fD*T_s; 
                        aa = 2*pi*fDTs*(T:-1:1);
                        bb = 2*pi*fDTs*(1:DataSym); 
                        for l = 1:L
                            %ACF for channel estimation
                            rho1 = (ricianFactor(k,l)/(ricianFactor(k,l)+1)) * exp(1i*aa*cos(Garma(k) - Theta(k,l)));
                            rho2 = (1/(ricianFactor(k,l)+1)) * besseli(0, sqrt( - aa.^2 + ka_UE^2 + 1i*2*aa*ka_UE*cos(Garma(k) - Theta(k,l)))) ./ besseli(0, ka_UE);  
                            TCC_Matrix_Pilot_t(k,l,:) = rho1.' + rho2.'; 
                            %ACF for data transmission
                            rho1 = (ricianFactor(k,l)/(ricianFactor(k,l)+1)) * exp(1i*bb*cos(Garma(k) - Theta(k,l)));
                            rho2 = (1/(ricianFactor(k,l)+1)) * besseli(0, sqrt( - bb.^2 + ka_UE^2 + 1i*2*bb*ka_UE*cos(Garma(k) - Theta(k,l)))) ./ besseli(0, ka_UE);  
                            TCC_Matrix_Data_d(k,l,:) = rho1.' + rho2.';   
                        end 
                    end  
                    % Design service AP indicator U matrix according to different handover policies  
                    if n_yipi == 0  
                        %Initialize U matrix at n_yipi 
                        [U_fuzzy1, set_index_delete_fuzzy1] = servingAPgainInit(K,L,channelGaindB,max_num_APserving);
                        [U_fuzzy2, set_index_delete_fuzzy2] = servingAPgainInit(K,L,channelGaindB,max_num_APserving);  
                        [U_LSFC, set_index_delete_LSFC] = servingAPgainInit(K,L,channelGaindB,max_num_APserving); 
                    else%using U matrix 
                        %Prepare to store the values of HOM
                        Th_Gain1 = zeros(K,1);  
                        for k=1:K  
                            if UE_flag(k) == 1 
                                %If at a road boundary, re-initialize U_fuzzy2 matrix
                                U_fuzzy1(k,:) = 0;  
                                [~, Index] = sort(channelGaindB(k,:), 'descend'); 
                                max_n_indices_new_fuzzy1 = Index(1:max_num_APserving); 
                                U_fuzzy1(k,max_n_indices_new_fuzzy1) = 1; 
                            else 
                                %If not at the road boundary, perform the HOM calculation
                                U_temp_old_fuzzy1 = find(U_fuzzy1(k,:)==1);
                                %for Algorithm 1
                                [~, min_APindex_old_fuzzy1] = min(channelGaindB(k,U_temp_old_fuzzy1));    
                                %Calculate the HOM using fuzzy logic
                                Th_Gain1(k,1)=evalfis(fis,[v_r(k),sum(abs(TCC_Matrix_Data_d(k,min_APindex_old_fuzzy1,:)))/DataSym,channelGaindB(k,min_APindex_old_fuzzy1)]);
                            end
                        end 

                        %U matrix and number of handover using fuzzy&Dir based method
                        [U_fuzzy1, set_index_delete_fuzzy1,Num_handover_kN_fuzzyBased_Dir(1:K,n_i,UE_lambda_i,serving_num_AP_i,n),Num_PP_handover_kN_fuzzyBased_Dir(1:K,n_i,UE_lambda_i,serving_num_AP_i,n)] = ...
                            fuzzyLogicDir(K,U_fuzzy1,set_index_delete_fuzzy1,max_num_APserving,channelGaindB,Th_Gain1,Garma,Theta);

                       %Prepare to store the values of HOM
                        Th_Gain2 = zeros(K,1);  
                        for k=1:K  
                            if UE_flag(k) == 1 
                                %If at a road boundary, re-initialize U_fuzzy2 matrix
                                U_fuzzy2(k,:) = 0;  
                                [~, Index] = sort(channelGaindB(k,:), 'descend'); 
                                max_n_indices_new_fuzzy2 = Index(1:max_num_APserving); 
                                U_fuzzy2(k,max_n_indices_new_fuzzy2) = 1; 
                            else
                                %If not at the road boundary, perform the HOM calculation
                                U_temp_old_fuzzy2 = find(U_fuzzy2(k,:)==1);
                                %for Algorithm 1
                                [~, min_APindex_old_fuzzy2] = min(channelGaindB(k,U_temp_old_fuzzy2));    
                                %Calculate the HOM using fuzzy logic
                                Th_Gain2(k,1)=evalfis(fis,[v_r(k),sum(abs(TCC_Matrix_Data_d(k,min_APindex_old_fuzzy2,:)))/DataSym,channelGaindB(k,min_APindex_old_fuzzy2)]);
                            end
                        end 

                        %U matrix and number of handover using fuzzy based method
                        [U_fuzzy2, set_index_delete_fuzzy2,Num_handover_kN_fuzzyBased_nonDir(1:K,n_i,UE_lambda_i,serving_num_AP_i,n),Num_PP_handover_kN_fuzzyBased_nonDir(1:K,n_i,UE_lambda_i,serving_num_AP_i,n)] = ...
                            fuzzyLogicNonDir(K,U_fuzzy2,set_index_delete_fuzzy2,max_num_APserving,channelGaindB,Th_Gain2);
                        %U matrix and number of handover using LSFC based method
                        [U_LSFC,set_index_delete_LSFC,Num_handover_kN_LSFC(1:K,n_i,UE_lambda_i,serving_num_AP_i,n),Num_PP_handover_kN_LSFC(1:K,n_i,UE_lambda_i,serving_num_AP_i,n)] = ...
                            maxGain(K,L,UE_flag,channelGaindB,max_num_APserving,U_LSFC,set_index_delete_LSFC);  
                    end   
                end
            end
        end 
    end
end
 
 
 