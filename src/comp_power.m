% Compute power for multiple normality tests: MK+WANG, MK+HZ_hS+WANG,
% MK+MS+HZ_hS+HZ_hL

clear
% clc
rng('shuffle')

n_samples=[20:20:100 200 400];
k_variate=[2 3 4 5 7 10];
N=50000; % Monte Carlo experiments
a = 0.05; % alpha
POWER_MW = zeros(length(n_samples), length(k_variate));
POWER_MBW = POWER_MW;
POWER_MMBB = POWER_MW;
POWER_ALL = zeros(length(n_samples), length(k_variate), 7);

for j=1:length(n_samples)
    n=n_samples(j);
    for p=1:length(k_variate)
        tic
        %--- load data ----
        fprintf('Processing for sample size n=%d...k=%d...\n',n,k_variate(p));
        str=strcat('new_data/Cnu_u_p_',num2str(k_variate(p)),'_n_',num2str(n));
        load(str); % variable Cnu
        str=strcat('new_data/Cnhu_p_',num2str(k_variate(p)),'_n_',num2str(n));
        load(str); % variable Cnhu
        u_MW=find(a < Cn_u(:,1), 1)-1;
        u_MBW=find(a < Cn_u(:,2), 1)-1;
        u_MMBB=find(a < Cn_u(:,3), 1)-1;
        CV_MW = Cnhu(u_MW,:);
        CV_MBW = Cnhu(u_MBW,:);
        CV_MMBB = Cnhu(u_MMBB,:);
        CV_all=Cnhu(a/0.0001,:);    %組合內個別檢定方法臨界值
        
        u1 = 0; u2 = 0; u3 = 0; u_all = zeros(1,7);
        for i=1:N
%-----------------------對立分配資料選擇   
%             X = multi_norm_sample_from(1,n_samples(j), k_variate(p), NaN); % Pearson II (m=0)
            X = multi_norm_sample_from(16,n_samples(j), k_variate(p)); % Pearson II (m=0)
    %      X=multi_T(p,0.5,n);
    %      X=multi_norm_sample_from(r,n,p,seed);
    %        X=multi_norm_sample_MVN(0.9,0,0,3,n,p,seed);
    %      X=multi_PeasonVII(p,10,0.5,n);
    %      X=multi_Khintchine1(p,n,seed);
%--------------------------------------------------
            [MK, MS, h141, hS, hL, WANG]=compute_ststs(X); % compute test statistic for each test
            Tmp = [MK MK WANG] - CV_MW([1 2 7]);
            cnt=max([Tmp(1).*Tmp(2) -Tmp(3) ], [], 2) > 0;
            u1 = u1 + cnt;
            
            Tmp = [MK MK hS WANG] - CV_MBW([1 2 5 7]);
            cnt=max([Tmp(1).*Tmp(2) Tmp(3) -Tmp(4)], [], 2) > 0;
            u2 = u2 + cnt;
            
            Tmp = [MK MK MS hS hL] - CV_MMBB([1 2 3 5 6]);
            cnt=max([Tmp(1).*Tmp(2) Tmp(3:end)], [], 2) > 0;
            u3 = u3 + cnt;
            
            Tmp = [MK MK MS h141 hS hL WANG] - CV_all;
            cnt=[Tmp(1).*Tmp(2) Tmp(3:end) -Tmp(end)] > [0 0 0 0 0 0 0];
            u_all = u_all + cnt;
        end
        fprintf('Time Elapsed %10.2f \n', toc);
        POWER_MW(j,p)=u1/N;
        POWER_MBW(j,p)=u2/N;
        POWER_MMBB(j,p)=u3/N;
        POWER_ALL(j,p,:)=u_all/N;
    end
end
% str=strcat('new_data/POWER_Normal_alpha_00',num2str(a*100));
str=strcat('new_data/POWER_BurrParetoLogistic(a1)_alpha_00',num2str(a*100));

save(str, 'POWER_MW', 'POWER_MBW', 'POWER_MMBB', 'POWER_ALL')

function [MK, MS, h141, hS, hL, WANG] = compute_ststs(X)
% This function computes and output 6 test statistics for a sample data X.

    MK = multi_norm_MK_CV(X, 1);
    MS = multi_norm_MS_CV(X, 1);
    [h141, hS, hL] = multi_norm_HZ_CV(X,1);
    WANG = multi_norm_percentile_simu(X,10000,0.05);
end
