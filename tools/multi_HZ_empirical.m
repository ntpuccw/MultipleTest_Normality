% Calculate critical values for HZ under various scenarios
% Output (stored in file HZ_h141_k_? ) contains three CV, which are
% HZ_h141, HZ_hS(BhS), HZ_hL(BhL)

clc
clear
rng('shuffle')

n_samples= [20:20:100 200 400];%input('Enter samples: ');%10:10:100;
k_variate=[2 3 4 5 7 10];%input('Enter k_varuate:' );%2:5;
N=10000;
% alfa=[0.001 0.01 0.025 0.05 0.075 0.1 0.2];
alfa = 0.0001:0.0001:0.3;
rep = [30 30 30 30 30 20 10]; % repeated times for each sample size
% rep_choices=[30 20 10 5];

for k=k_variate
    HZ_h141=zeros(length(n_samples),length(alfa)+1);
    HZ_hS=zeros(length(n_samples),length(alfa)+1);
    HZ_hL=zeros(length(n_samples),length(alfa)+1);
    cnt=1;
    for j=1:length(n_samples)
%         rep=rep_choices(ceil(n/50));
        HZ_h141_tmp=zeros(rep(j),length(alfa));
        HZ_hS_tmp=zeros(rep(j),length(alfa));
        HZ_hL_tmp=zeros(rep(j),length(alfa));
        tic
        for i=1:rep(j)
%             Y=multi_norm_sample_from(0,n,N,k);
%             HZ_h141_stats=zeros(N,1);
            HZ_stats=zeros(N,3);
        %----------------------
            fprintf('Rep=%3.0f    Sample size  n=%3.0f ... variates k=%2.0f \n',i,n_samples(j),k)
            for s=1:N
                X = multi_norm_sample_from(0,n_samples(j),k,NaN);
%                 X=Y(((i-1)*n+1):i*n,:);%n x k
                [HZ_stats(s,1), HZ_stats(s,2), HZ_stats(s,3)]=multi_norm_HZ_CV(X,1);
            end       
                [YCDF,XCDF] = ecdf(HZ_stats(:,1));
                HZ_h141_tmp(i,:)= XCDF(round(length(XCDF)*(1-alfa)))';
                [YCDF,XCDF] = ecdf(HZ_stats(:,2));
                HZ_hS_tmp(i,:)= XCDF(round(length(XCDF)*(1-alfa)))';
                [YCDF,XCDF] = ecdf(HZ_stats(:,3));
                HZ_hL_tmp(i,:)= XCDF(round(length(XCDF)*(1-alfa)))';
        end
        HZ_h141(cnt,:)=[n_samples(j) mean(HZ_h141_tmp)];
        HZ_hS(cnt,:)=[n_samples(j) mean(HZ_hS_tmp)];
        HZ_hL(cnt,:)=[n_samples(j) mean(HZ_hL_tmp)];
        cnt=cnt+1;
        t=toc;
        fprintf('Time Elapsed %10.2f \n', t);
    end
    str=strcat('save new_data/HZ_h141_k',num2str(k), ' HZ_h141 HZ_hS HZ_hL  alfa');
    eval(str);
end
