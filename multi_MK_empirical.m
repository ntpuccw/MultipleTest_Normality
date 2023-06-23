% Calculate critical values for Madias's Multivariate Kurtosis b2p under various scenarios
% The outputs are stored in file MK_k_? for various k dimensions.
% The output variable MK_LOWER, MK_UPPER represent the lower CV and upper
% CV for the b2p statistic (a two-tailed test)

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
    MK_UPPER=zeros(length(n_samples),length(alfa)+1);
    MK_LOWER=zeros(length(n_samples),length(alfa)+1);
    cnt=1;
    for j=1:length(n_samples)
%         rep=rep_choices(ceil(n/50));
        MK_tmp_upper=zeros(rep(j),length(alfa));
        MK_tmp_lower=zeros(rep(j),length(alfa));
        tic
        for i=1:rep(j)
%             Y=multi_norm_sample_from(0,n,N*k,sum(100*clock));
            MK_b2p=zeros(N,1);
        %----------------------
            fprintf('Rep=%d    Sample size  n=%d ... variates k=%d \n',i,n_samples(j),k)
            for s=1:N
                X = multi_norm_sample_from(0,n_samples(j),k,NaN);
%                 X=Y(((i-1)*n+1):i*n,:);%n x k
                MK_b2p(s)=multi_norm_MK_CV(X,1);
            end       
            [YCDF,XCDF] = ecdf(MK_b2p);
            MK_tmp_upper(i,:)= XCDF(round(length(XCDF)*(1-alfa/2)))';
            pos=round(length(XCDF)*alfa/2);
            pos(pos==0)=1;
            MK_tmp_lower(i,:)= XCDF(pos)';
        end
        MK_UPPER(cnt,:)=[n_samples(j) mean(MK_tmp_upper)];
        MK_LOWER(cnt,:)=[n_samples(j) mean(MK_tmp_lower)];
        cnt=cnt+1;
        t=toc;
        fprintf('Time Elapsed %10.2f \n', t);
    end
    str=strcat('save new_data/MK_k',num2str(k), ' MK_UPPER MK_LOWER alfa');
    eval(str);
end
