%The script can be used to show the number of handover 
%according to the soft handover algorithm in the manuscript:
 
%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
% None
  
%Output:
% Show the results of <numHandoverDiffAPnum.m> 
 
clc 

Num_handover_kN_fuzzyBased_Dir_result = squeeze(mean(mean(sum(Num_handover_kN_fuzzyBased_Dir,1),2),4));
Num_PP_handover_kN_fuzzyBased_Dir_result = squeeze(mean(mean(sum(Num_PP_handover_kN_fuzzyBased_Dir,1),2),4));
Num_handover_kN_fuzzyBased_nonDir_result = squeeze(mean(mean(sum(Num_handover_kN_fuzzyBased_nonDir,1),2),4));
Num_PP_handover_kN_fuzzyBased_nonDir_result = squeeze(mean(mean(sum(Num_PP_handover_kN_fuzzyBased_nonDir,1),2),4));
Num_handover_kN_gainBased_result = squeeze(mean(mean(sum(Num_handover_kN_LSFC,1),2),4));
Num_PP_handover_kN_gainBased_result = squeeze(mean(mean(sum(Num_PP_handover_kN_LSFC,1),2),4)); 
 
for num_AP_i = 1:length(num_AP_number) 
    if num_AP_i == 1
        data1 = Num_handover_kN_gainBased_result(num_AP_i);
        data2 = Num_PP_handover_kN_gainBased_result(num_AP_i); 
        data3 = Num_handover_kN_fuzzyBased_nonDir_result(num_AP_i);
        data4 = Num_PP_handover_kN_fuzzyBased_nonDir_result(num_AP_i); 
        data5 = Num_handover_kN_fuzzyBased_Dir_result(num_AP_i);
        data6 = Num_PP_handover_kN_fuzzyBased_Dir_result(num_AP_i);
    else 
        data1 = [data1 Num_handover_kN_gainBased_result(num_AP_i)];
        data2 = [data2 Num_PP_handover_kN_gainBased_result(num_AP_i)];  
        data3 = [data3 Num_handover_kN_fuzzyBased_nonDir_result(num_AP_i)];
        data4 = [data4 Num_PP_handover_kN_fuzzyBased_nonDir_result(num_AP_i)]; 
        data5 = [data5 Num_handover_kN_fuzzyBased_Dir_result(num_AP_i)];
        data6 = [data6 Num_PP_handover_kN_fuzzyBased_Dir_result(num_AP_i)]; 
    end 
end
 
figure
box on, grid on, hold on
subplot(2,2,1)
grid on, box on, hold on
b = bar([num_AP_number],[data1; data3; data5], 'grouped'); 
ylim([0 50])
xlabel('$L$','Interpreter','Latex')
ylabel('$\mathrm{N_{H}}$','Interpreter','Latex')
legend( 'LSFC based','Proposed alg. (fuzzy)','Proposed alg. (fuzzy&dir.)')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';
ylim([0 80])

subplot(2,2,2)
grid on, box on, hold on
b = bar([num_AP_number],[data2; data4; data6], 'grouped');   
xlabel('$L$','Interpreter','Latex')
ylabel('$\mathrm{N_{PPH}}$','Interpreter','Latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';
ylim([0 20])

subplot(2,2,[3,4])
grid on, box on, hold on
b = bar([num_AP_number],[data2./data1; data4./data3; data6./data5].*100, 'grouped');  
xlabel('$L$','Interpreter','Latex')
ylabel('PPH Radio (%)','Interpreter','Latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';
ylim([0 45])
 