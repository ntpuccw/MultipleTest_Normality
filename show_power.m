clear
dirpath = "new_data/";
filename =[ "POWER_T(10)_alpha_005";
            "POWER_KHN(normal)_alpha_005";
            "POWER_KHN(exp)_alpha_005";
            "POWER_PearsonVII(m10)_alpha_005";
            "POWER_PearsonII(m10)_alpha_005";
            "POWER_PearsonII(m05)_alpha_005";
            "POWER_PearsonII(m0)_alpha_005";
            "POWER_SkewNormal(a2)_alpha_005";
            "POWER_SkewNormal(a4)_alpha_005";
            "POWER_MixedMVN(05B)_alpha_005";
            "POWER_MixedMVN(05mu3)_alpha_005";
            "POWER_MixedNormal(09B)_alpha_005";
            "POWER_MixedNormal(09mu3)_alpha_005";
            "POWER_MixedNormal(079mu3)_alpha_005";
            "POWER_Normal_alpha_001";
            "POWER_Normal_alpha_005"];
fIdx = 14;        
load(dirpath+filename(fIdx))
n = [20:20:100 200 400];
d = [2 3 4 5 7 10];

for ni=1:length(n)
% ni = 5;
    figure
    plot(POWER_MW(ni,:),'-*', 'DisplayName', 'MW')
    hold on
    plot(POWER_MBW(ni,:),'--^', 'DisplayName', 'MBW')
    plot(POWER_MMBB(ni,:),':s', 'DisplayName', 'MMBB')
%     plot(POWER_ALL(ni,:,1),'-.o', 'DisplayName', 'MK') % MK
    % plot(POWER_ALL(ni,:,2),'-.+', 'DisplayName', 'MS') % MS
    % plot(POWER_ALL(ni,:,3),'-.|', 'DisplayName', 'HZ') % H141
    % plot(POWER_ALL(ni,:,4),'-.x', 'DisplayName', 'HS') % HS
    % plot(POWER_ALL(ni,:,5),'-.<', 'DisplayName', 'HL') % HL
    % plot(POWER_ALL(ni,:,7),'-.>', 'DisplayName', 'WANG') % WANG
    hold off
    title("n = " + num2str(n(ni)))
    xticklabels(d)
    ylabel('empirical power')
    xlabel('dimension')
    grid on
    legend('Location','northwest')
    if fIdx == 15 
        ylim([0.00 0.02])
    elseif fIdx == 16 
        ylim([0.04 0.06])
    else
        ylim([0 1])
    end
end