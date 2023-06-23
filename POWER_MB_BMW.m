clc
clear
% �Ω�p��X�h���˩w�զX MB �M BMW ���˩w�O�A��i�H�p��զX���C�ؤ�k�˩w�O
% MB  �h���˩w�զX�]�t: MS�BMK�BB(h_S)�BB(h_L)   Tenreiro (2011)
% BMW �h���˩w�զX�]�t:MK�BB(h_S)�BW_{min}(5)    Wang and Hwang (2011)
% �u���զX���@�ӥH�W�ڵ���L���]�A�h�h���˩w�ڵ���L���]
% MS�BB(h_S)�BB(h_L) �k���˩w�A�˩w�έp�q-���{�ɭ�>0 �ڵ���L���]
% MK �����˩w�A(MK-�{�ɭȤW��)*(MK-�{�ɭȤU��)>0 �ڵ���L���]
% W_{min}(5) �����˩w�A"�t"�˩w�έp�q-"�t"�{�ɭ�>0 �ڵ���L���]�A�GCV�Ȭҥ[�W"�t"
%
% �{���c�y:
%   1.�q MBtest_u �����ӧO���@�~�t u �ϱo���髬�@�~<=alpha �A�ñq T �����������{�ɭ�
%   2.�i�� N �չ�ߤ��t�h�ܶq�`�A�˩w MB �M BMW 
%   3.�˩w���G outputs :POWER_MB �M POWER_BMW  (�a�b n_sample, ��bk_variate)
%   4.�@��
% 
%created by �x�_�j�ǲέp�Ǩt�@��s�͡G���P��  2014.7
%---------------------�]�w
n_sample=20;  % �{�ɭȸ�Ʀ� n=20,60,100,400
k_variate=[2 3 4 5 10];% �{�ɭȸ�Ʀ� d=2,3,4,5,10
alpha=0.05;
N=50000;
seed=781103;
POWER_BMW=zeros(4,10);
POWER_MB=zeros(4,10);
% POWER=zeros(20,10);
r=1;                  % multi_norm_sample_from ��ߤ��t���

for j=1:length(n_sample)
    n=n_sample(j);
    for p=k_variate
        tic
        %-------------------------Ū�����@�~�t' MBtest_u'�Ωҹ������{�ɭ�' T'
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
        %CV3=T(alpha/0.0001,1:6); %�զX���ӧO�˩w��k�{�ɭ�
        
        %----------------�p�ƾ��k0
        u1=0;
        u2=0;
%       u=zeros(1,5);
        
        for i=1:N
%-----------------------��ߤ��t��ƿ��     
%      X=multi_T(p,0.5,n);
%      X=multi_norm_sample_from(r,n,p,seed);
       X=multi_norm_sample_MVN(0.9,0,0,3,n,p,seed);
%      X=multi_PeasonVII(p,10,0.5,n);
%      X=multi_Khintchine1(p,n,seed);
%--------------------------------------------------
     
            X=(chol(cov(X))'\X')';
            [b2p,MS,HS,HL] = combi_STAT(X);
            E=-multi_norm_percentile_simu(X,10000,0.05);     %�����˩w�A���F��K �˩w�έp�q �M �{�ɭ� �ҥ[�W"�t"��
            
            T=[b2p b2p MS HS HL E]-CV_MB;
            cnt1=max([T(1)*T(2) T(3) T(4) T(5)])>0;
            u1=u1+cnt1;% MB �զX�˩w�O
            
            T=[b2p b2p MS HS HL E]-CV_BMW;
            cnt2=max([T(1)*T(2) T(4) T(6)])>0;
            u2=u2+cnt2;% BMW �զX�˩w�O
            
%             T=[b2p b2p MS HS HL E]-CV3;
%             cnt=([T(1)*T(2) T(3) T(4) T(5) T(6)])>[0 0 0 0 0];
%             u=u+cnt;% �ӧO�˩w��k�˩w�O
            
            seed=seed+1;

        end
        
        t=toc;
        fprintf('Time Elapsed %10.2f \n', t);
        POWER_MB(j,p)=u1/N;
        POWER_BMW(j,p)=u2/N;
%       POWER(5*j-4:5*j,p)=u'./N;
    end
end




%--------------------------------------------------�@��
n=20;  % ���e�˥��Ƭ��h��
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
