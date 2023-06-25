% 計算論文表 2、表 3 的 U_n,a 值
clear

n_samples=[20 60 100 200 400];
k_variate=[2 3 4 5 7 10];

a = 0.05; % alpha
u_na_MW = zeros(length(n_samples), length(k_variate));
u_na_MBW = u_na_MW;
u_na_MMBB = u_na_MW;


for j=1:length(n_samples)
    n=n_samples(j);
    for p=1:length(k_variate)
        %--- load data ----
        str=strcat('new_data/Cnu_u_p_',num2str(k_variate(p)),'_n_',num2str(n));
        load(str); % variable Cnu
        u_na_MW(j,p) = find(a < Cn_u(:,1), 1)-1;
        u_na_MBW(j,p) = find(a < Cn_u(:,2), 1)-1;
        u_na_MMBB(j,p) = find(a < Cn_u(:,3), 1)-1;
    end
end
table_MW = u_na_MW*0.0001
table_MBW = u_na_MBW*0.0001
table_MMBB = u_na_MMBB*0.0001