function Y=multi_norm_sample_from(option,n,p,seed)
% This function generates multivariate data with prescribed dist.
%
% syntax : Y=multi_norm_sample_from(option,n,p)
% inputs :
%   option: indicate which distribution to generate from
%   n: sample size
%   p: p-variate
%   seed:  set the random state
%
% outputs :
%   Y: multivariate data
%
% created by    Chun-Chao Wang
%               Dept. of Statistics    
%               National Taipei University
%               Taipei, Taiwan
%
% References :
% Wang CC, Hwang YT (2010). "A new functional statistic for multivariate normality." 
%       Statistics and Computing, 21(4), 501-509..
% =================================================================
% if nargin<4
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% else
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed', seed));
% end
if isnan(seed) % user does not enter any seed, use a random seed
    myseed=sum(100*clock);
else % user seed to repeat the same computation
    myseed=seed;
end
if verLessThan('matlab','8.0')
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',myseed));
else
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',myseed));
end    

N=1;% reserved for number of Monte Carlo replications.
    switch option
        case 0 % H0: Multivariate Normal
            rho=0.5;% correlation
            mu=zeros(1,p);
            Sigma=eye(p);Sigma(Sigma==0)=rho;%equal rho
            R=chol(Sigma);
            Y=repmat(mu,n,1) + randn(n,p)*R;
%             Y=mvnrnd(mu, Sigma, n*N);          % needs statistics toolbox  
        case 1  %multivariate Uniform
            m=0;% Pearson Type II parameter
            Y=multi_PearsonII(p,m,n*N);            
        case 2 %multivariate T(3)
            nu=3;% Multivariate T parameter
            rho=0;
            Y=multi_T(p,nu,rho,n*N);            
        case 3 %multivariate T(10)
            nu=10;% Multivariate T parameter
            rho=0;
            Y=multi_T(p,nu,rho,n*N);            

        case 4 % multivariate Khintchine
%             rho=0;      %the correlation among  Ui
            Y=multi_Khintchine(p,n*N);
                            
        case 5 %mixed multivariate normal
            pi1=0.1;
%             mu1=zeros(1,p);
%             Sigma1=eye(p);%equal rho            
            mu2=0*ones(1,p);
            Sigma2=eye(p);
            Sigma2(Sigma2==0)=0.9;% off diagonal
            Sigma2(Sigma2==1)=1;%main diagonal
            N1=binornd(n,pi1,1,N);
            N2=n-N1;
            Y1= randn(sum(N1),p);% standard normal
            R2=chol(Sigma2);
            Y2=repmat(mu2,sum(N2),1) + randn(sum(N2),p)*R2;            
%             Y1=mvnrnd(mu1, Sigma1, sum(N1));
%             Y2=mvnrnd(mu2, Sigma2, sum(N2)); 
            Y=zeros(n*N,p);
            N1=[0 cumsum(N1)];
            N2=[0 cumsum(N2)];
            for i=1:N
                Y((i-1)*n+1:i*n,:)=[Y1(N1(i)+1:N1(i+1),:);Y2(N2(i)+1:N2(i+1),:)];
            end
            
        case 6%multivariate skew nomal with mild skewness
            a=2*ones(1,p);
            rho=0.5;
            mu=zeros(1,p);
            Sigma=eye(p);Sigma(Sigma==0)=rho;
            Y = rmsn(n*N,mu,Sigma,a);
        case 7 %multivariate skew nomal with moderate skewness
            a=4*ones(1,p);
            rho=0.5;
            mu=zeros(1,p);
            Sigma=eye(p);Sigma(Sigma==0)=rho;
            Y = rmsn(n*N,mu,Sigma,a);
    end