clc
clear
% 用於計算出多重檢定組合 MB 和 BMW 的檢定力，亦可以計算組合內每種方法檢定力
% MB  多重檢定組合包含: MS、MK、B(h_S)、B(h_L)   Tenreiro (2011)
% BMW 多重檢定組合包含:MK、B(h_S)、W_{min}(5)    Wang and Hwang (2011)
% 只有組合內一個以上拒絕虛無假設，則多重檢定拒絕虛無假設
% MS、B(h_S)、B(h_L) 右尾檢定，檢定統計量-該臨界值>0 拒絕虛無假設
% MK 雙尾檢定，(MK-臨界值上界)*(MK-臨界值下界)>0 拒絕虛無假設
% W_{min}(5) 左尾檢定，"負"檢定統計量-"負"臨界值>0 拒絕虛無假設，故CV值皆加上"負"
%
% 程式構造:
%   1.從 MBtest_u 中找到個別型一誤差 u 使得整體型一誤<=alpha ，並從 T 中找到對應的臨界值
%   2.進行 N 組對立分配多變量常態檢定 MB 和 BMW 
%   3.檢定結果 outputs :POWER_MB 和 POWER_BMW  (縱軸 n_sample, 行軸k_variate)
%   4.作圖
% 
%created by 台北大學統計學系　研究生：陳燊飛  2014.7
%---------------------設定
n_sample=20;  % 臨界值資料有 n=20,60,100,400
k_variate=[2 3 4 5 10];% 臨界值資料有 d=2,3,4,5,10
alpha=0.05;
N=50000;
seed=781103;
POWER_BMW=zeros(4,10);
POWER_MB=zeros(4,10);
% POWER=zeros(20,10);
r=1;                  % multi_norm_sample_from 對立分配選擇

for j=1:length(n_sample)
    n=n_sample(j);
    for p=k_variate
        tic
        %-------------------------讀取型一誤差' MBtest_u'及所對應的臨界值' T'
        fprintf('Processing for sample size n=%3.0f..k=%3.0f...%3.0f...\n',n,p);
        str=strcat('load data/MBtest_u_u_p_',num2str(p),'_n_',num2str(n), ' MBtest_u');
        eval(str);
        str=strcat('load data/MBtest_CV_CV_p_',num2str(p),'_n_',num2str(n), ' T');
        eval(str);
        W_MB=find(alpha<MBtest_u(:,1));
        W_NEW=find(alpha<MBtest_u(:,2));
        W1=W_MB(1)-1;
        W2=W_NEW(1)-1;
        CV_MB=T(W1,:);
        CV_BMW=T(W2,:);
        %CV3=T(alpha/0.0001,1:6); %組合內個別檢定方法臨界值
        
        %----------------計數器歸0
        u1=0;
        u2=0;
%       u=zeros(1,5);
        
        for i=1:N
%-----------------------對立分配資料選擇     
%      X=multi_T(p,0.5,n);
%      X=multi_norm_sample_from(r,n,p,seed);
       X=multi_norm_sample_MVN(0.9,0,0,3,n,p,seed);
%      X=multi_PeasonVII(p,10,0.5,n);
%      X=multi_Khintchine1(p,n,seed);
%--------------------------------------------------
     
            X=(chol(cov(X))'\X')';
            [b2p,MS,HS,HL] = combi_STAT(X);
            E=-multi_norm_percentile_simu(X,10000,0.05);     %左尾檢定，為了方便 檢定統計量 和 臨界值 皆加上"負"號
            
            T=[b2p b2p MS HS HL E]-CV_MB;
            cnt1=max([T(1)*T(2) T(3) T(4) T(5)])>0;
            u1=u1+cnt1;% MB 組合檢定力
            
            T=[b2p b2p MS HS HL E]-CV_BMW;
            cnt2=max([T(1)*T(2) T(4) T(6)])>0;
            u2=u2+cnt2;% BMW 組合檢定力
            
%             T=[b2p b2p MS HS HL E]-CV3;
%             cnt=([T(1)*T(2) T(3) T(4) T(5) T(6)])>[0 0 0 0 0];
%             u=u+cnt;% 個別檢定方法檢定力
            
            seed=seed+1;

        end
        
        t=toc;
        fprintf('Time Elapsed %10.2f \n', t);
        POWER_MB(j,p)=u1/N;
        POWER_BMW(j,p)=u2/N;
%       POWER(5*j-4:5*j,p)=u'./N;
    end
end




%--------------------------------------------------作圖
n=20;  % 欲畫樣本數為多少
j=find(n_sample==n);
X=0:4;
Y1=POWER_MB(j,[2 3 4 5 10]);
plot(X,Y1,'-.*k')
hold on
Y2=POWER_BMW(j,[2 3 4 5 10]);
plot(X,Y2,'-*k')

% Y=POWER(5*j-4:5*j,[2 3 4 5 10]);
% plot(X,Y(1,:),'-.sk')
% plot(X,Y(2,:),'-sk')
% plot(X,Y(3,:),'-^k')
% plot(X,Y(4,:),'-.^k')
% plot(X,Y(5,:))
% h = legend('MB','BMW','MK','MS','B(h_S)','B(h_L)','W_{min,m}(5)');
h = legend('MB','BMW');
set(h,'Interpreter','none')
set(gca,'XTickLabel',{'2','','3','','4','','5','','10'})
ylim([0 1])
xlabel('d')
ylabel('power')
