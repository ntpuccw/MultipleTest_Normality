clear
rng(0)
n=10;
p=2;
pi1=0.1;
N=1;
%             mu1=zeros(1,p);
%             Sigma1=eye(p);%equal rho     
NN=1
tic
for j=1:NN
            mu1=0*ones(1,p);
            Sigma1=eye(p);
            mu2=0*ones(1,p);
            Sigma2=eye(p);
            Sigma2(Sigma2==0)=0.9;% off diagonal
            N1=binornd(n,pi1);
            N2=n-N1;
            Y1=repmat(mu1,N1,1) + randn(N1,p);% standard normal
            R2=chol(Sigma2);
            Y2=repmat(mu2,N2,1) + randn(N2,p)*R2;            
            Y = [Y1;Y2];
%             Y=zeros(n*N,p);
%             N1=[0 cumsum(N1)];
%             N2=[0 cumsum(N2)];
%             for i=1:N
%                 Y((i-1)*n+1:i*n,:)=[Y1(N1(i)+1:N1(i+1),:);Y2(N2(i)+1:N2(i+1),:)];
%             end
end         
   toc
   Y

   tic
   for j=1:NN
            rng(0)
            mu1=zeros(1,p);
            Sigma1=eye(p);%equal rho    
            N1=binornd(n,pi1,1,N);
            N2=n-N1;
            Y1=mvnrnd(mu1, Sigma1, sum(N1));
            Y2=mvnrnd(mu2, Sigma2, sum(N2)); 
            YY = [Y1;Y2]
   end
toc
YY
rng(0)
YYY=mixed_MVN(mu1, mu2, Sigma2, n, p, pi1) 


function Y = mixed_MVN(mu1, mu2, Sigma2, n, p, pi1)    
        N1=binornd(n,pi1);
        N2=n-N1;
        Y1=repmat(mu1,N1,1) + randn(N1,p);% standard normal
        R2=chol(Sigma2);
        Y2=repmat(mu2,N2,1) + randn(N2,p)*R2;            
        Y = [Y1;Y2];
    end
            