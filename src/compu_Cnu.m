% Compute Cn,u for each u = 0.0001:0.0001:0.3
% The program will load pre-stored data Cnhu_p_?_n_? for p=2,3,4,5,7,10, n =[20:20:100 200 400]
% In each Cnhu_p_?_n_?, a variable Cnhu is used, where 
% Cnhu=[MK_LOWER MK_UPPER MS HZ141 HZ_hS HZ_hL WANG] stores CV for each
% single test.
% This program compute the type I errors for three combination of single tests
% and stores the result variable Cn_u in Cnu_u_p_?_n_?.mat.
% where Cn_u = [MK+WANG MK+HZ_hS+WANG MK+MS+HZ_hs+HZ_hL]
% each value represents the type I error in 50000 simulations.
clear
% clc
rng('shuffle')

n_samples= [20:20:100 200 400];
k_variate=[2 3 4 5 7 10];
N=50000; % Monte Carlo experiments

for p=k_variate
    for j=1:length(n_samples)
        tic
        fprintf('Processing for sample size n=%d..k=%d ...\n',n_samples(j),p);
        str=strcat('new_data/Cnhu_p_',num2str(p),'_n_',num2str(n_samples(j)));
        load(str); % variable Cnhu
        Cnhu_MW=Cnhu(:,[1 2 7]); % MK + WANG
        Cnhu_MBW=Cnhu(:,[1 2 5 7]); % MK + BhS + WANG
        Cnhu_MMBB=Cnhu(:,[1 2 3 5 6]); % MK + MS + BhS + BhL
        
        nu = size(Cnhu,1);
        u1 = zeros(nu,1);
        u2 = zeros(nu,1);
        u3 = zeros(nu,1);
          %----------------
        for i=1:N
            X = multi_norm_sample_from(0,n_samples(j), p, NaN);
            [MK, MS, h141, hS, hL, WANG]=compute_ststs(X); % compute test statistic for each test
            CV_MW = repmat([MK MK WANG],nu,1);
            CV_MBW = repmat([MK MK hS WANG],nu,1);
            CV_MMBB = repmat([MK MK MS hS hL],nu,1);
            
            Tmp = CV_MW - Cnhu_MW;
            tmp1=max([Tmp(:,1).*Tmp(:,2) -Tmp(:,3) ], [], 2);
            Tmp=CV_MBW - Cnhu_MBW;
            tmp2=max([Tmp(:,1).*Tmp(:,2) Tmp(:,3) -Tmp(:,4)], [], 2);
            Tmp=CV_MMBB - Cnhu_MMBB;
            tmp3=max([Tmp(:,1).*Tmp(:,2) Tmp(:,3:end)], [], 2);
            cnt1 = tmp1 > 0;
            cnt2 = tmp2 > 0;
            cnt3 = tmp3 > 0;
            u1 = u1 + cnt1;
            u2 = u2 + cnt2;
            u3 = u3 + cnt3;
        end
        Cn_u=[u1 u2 u3]/N;
        str=strcat('new_data/Cnu_u_p_',num2str(p),'_n_',num2str(n_samples(j)));
        save(str, 'Cn_u')
        fprintf('Time Elapsed %10.2f \n', toc);   
    end
end

function [MK, MS, h141, hS, hL, WANG] = compute_ststs(X)
% This function computes and output 6 test statistics for a sample data X.

    MK = multi_norm_MK_CV(X, 1);
    MS = multi_norm_MS_CV(X, 1);
    [h141, hS, hL] = multi_norm_HZ_CV(X,1);
    WANG = multi_norm_percentile_simu(X,10000,0.05);
end