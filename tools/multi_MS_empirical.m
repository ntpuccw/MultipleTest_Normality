% % Calculate critical values for Madias's Multivariate Skewness b1p under various scenarios
% The outputs are stored in file MS_k_? for various k dimensions.
% The output variable MS represent CV for the b1p statistic (a upper-tailed test)
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
    MS=zeros(length(n_samples),length(alfa)+1);
    cnt=1;
    for j=1:length(n_samples)
%         rep=rep_choices(ceil(n/50));
        MS_tmp=zeros(rep(j),length(alfa));
        tic
        for i=1:rep(j)
%             Y=multi_norm_sample_from(0,n,N,k);
            MS_b1p=zeros(N,1);
        %----------------------
            fprintf('Rep=%3.0f    Sample size  n=%3.0f ... variates k=%2.0f \n',i,n_samples(j),k)
            for s=1:N
                X = multi_norm_sample_from(0,n_samples(j),k,NaN);
%                 X=Y(((i-1)*n+1):i*n,:);%n x k
                MS_b1p(s)=multi_norm_MS_CV(X,1);
            end       
                [YCDF,XCDF] = ecdf(MS_b1p);
                MS_tmp(i,:)= XCDF(round(length(XCDF)*(1-alfa)))';
        end
        MS(cnt,:)=[n_samples(j) mean(MS_tmp)];
        cnt=cnt+1;
        t=toc;
        fprintf('Time Elapsed %10.2f \n', t);
    end
    str=strcat('save new_data/MS_k',num2str(k), ' MS alfa');
    eval(str);
end
