%The script can be used to show the number of handover 
%according to the soft handover algorithm in the manuscript:
 
%Huafu Li, Yang Wang, Chenyang Sun, and Zhenyong Wang, "User-Centric
%Cell-Free Massive MIMO for IoT in Highly Dynamic Environments", submitted
%to IoTJ on May 29th, 2023.

%Input:
% None
  
%Output:
% Show the results of <numHandoverDiffUserDensity.m>
clc

Num_handover_kN_fuzzyBased_Dir_result = squeeze(mean(mean(sum(Num_handover_kN_fuzzyBased_Dir,1),2),5));
Num_PP_handover_kN_fuzzyBased_Dir_result = squeeze(mean(mean(sum(Num_PP_handover_kN_fuzzyBased_Dir,1),2),5));
Num_handover_kN_fuzzyBased_nonDir_result = squeeze(mean(mean(sum(Num_handover_kN_fuzzyBased_nonDir,1),2),5));
Num_PP_handover_kN_fuzzyBased_nonDir_result = squeeze(mean(mean(sum(Num_PP_handover_kN_fuzzyBased_nonDir,1),2),5));
Num_handover_kN_gainBased_result = squeeze(mean(mean(sum(Num_handover_kN_LSFC,1),2),5));
Num_PP_handover_kN_gainBased_result = squeeze(mean(mean(sum(Num_PP_handover_kN_LSFC,1),2),5));
   
for num_i = 1:length(max_num_APserving_range) 
    if num_i == 1
        data1 = Num_handover_kN_gainBased_result(num_i);
        data2 = Num_PP_handover_kN_gainBased_result(num_i); 
        data3 = Num_handover_kN_fuzzyBased_nonDir_result(num_i);
        data4 = Num_PP_handover_kN_fuzzyBased_nonDir_result(num_i); 
        data5 = Num_handover_kN_fuzzyBased_Dir_result(num_i);
        data6 = Num_PP_handover_kN_fuzzyBased_Dir_result(num_i); 
    else 
        data1 = [data1 Num_handover_kN_gainBased_result(num_i)];
        data2 = [data2 Num_PP_handover_kN_gainBased_result(num_i)];  
        data3 = [data3 Num_handover_kN_fuzzyBased_nonDir_result(num_i)];
        data4 = [data4 Num_PP_handover_kN_fuzzyBased_nonDir_result(num_i)]; 
        data5 = [data5 Num_handover_kN_fuzzyBased_Dir_result(num_i)];
        data6 = [data6 Num_PP_handover_kN_fuzzyBased_Dir_result(num_i)]; 
    end 
end

figure
box on, grid on, hold on 
subplot(2,2,1)
grid on, box on, hold on
bar([max_num_APserving_range],[data1; data3; data5], 'grouped')
xlabel('$L_s$','Interpreter','Latex')
ylabel('$\mathrm{N_{H}}$','Interpreter','Latex')
legend( 'LSFC based','Proposed (fuzzy)','Proposed (fuzzy&dir.)')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';  
ylim([0 100])

subplot(2,2,2)
grid on, box on, hold on
bar([max_num_APserving_range],[data2; data4; data6], 'grouped')   
xlabel('$L_s$','Interpreter','Latex')
ylabel('$\mathrm{N_{PPH}}$','Interpreter','Latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman'; 
ylim([0 30])

subplot(2,2,[3,4])
grid on, box on, hold on
bar([max_num_APserving_range],[data2./data1; data4./data3; data6./data5].*100, 'grouped')   
xlabel('$L_s$','Interpreter','Latex')
ylabel('PPH Radio (%)','Interpreter','Latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman'; 
ylim([0 45])



 
