% % Calculate critical values for Wang's Wmin(5) statistics under various scenarios
% The outputs are stored in file WANG_k_? for various k dimensions.
% The output variable WANG represent the  CV for WANG's Wmin(5) statistic 
% which is a lower-tailed test.

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
m = 10000;
q = 0.05;
for k=k_variate
    WANG=zeros(length(n_samples),length(alfa)+1);
    cnt=1;
    for j=1:length(n_samples)
%         rep=rep_choices(ceil(n/50));
        WANG_tmp=zeros(rep(j),length(alfa));
        tic
        for i=1:rep(j)
%             Y=multi_norm_sample_from(0,n,N,k);
            WANG_Wmin=zeros(N,1);
        %----------------------
            fprintf('Rep=%3.0f    Sample size  n=%3.0f ... variates k=%2.0f \n',i,n_samples(j),k)
            for s=1:N
                X = multi_norm_sample_from(0,n_samples(j),k,NaN);
%                 X=Y(((i-1)*n+1):i*n,:);%n x k
                WANG_Wmin(s)=multi_norm_percentile_simu(X,m,q);
            end       
                [YCDF,XCDF] = ecdf(WANG_Wmin);
                WANG_tmp(i,:)= XCDF(round(length(XCDF)*alfa))';
        end
        WANG(cnt,:)=[n_samples(j) mean(WANG_tmp)];
        cnt=cnt+1;
        t=toc;
        fprintf('Time Elapsed %10.2f \n', t);
    end
    str=strcat('save new_data/WANG_k',num2str(k), ' WANG alfa');
    eval(str);
end
