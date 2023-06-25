% Prepare C_n,h(u) from MK_lower, MK_upper, MS, BhS, BhL, Wmin(5)
% for various scenarios
% Each value in the matrix Cnhu represents the CV for each single test (column) 
% and various typr I error(row) from u=0.0001:0.0001:0.3

clear
n_samples= [20:20:100 200 400];%input('Enter samples: ');%10:10:100;
k_variate=[2 3 4 5 7 10];%input('Enter k_varuate:' );%2:5;


for p = k_variate
    str=strcat('new_data/MK_k',num2str(p));
    load(str)
    str=strcat('new_data/MS_k',num2str(p));
    load(str)
    str=strcat('new_data/HZ_h141_k',num2str(p));
    load(str)
    str=strcat('new_data/WANG_k',num2str(p));
    load(str)   
    nu = length(alfa); % number of u: 0.0001:0.0001:0.3
    Cnhu=zeros(nu, 7);
    for i = 1:length(n_samples)
        row = find(MK_LOWER(:,1)==n_samples(i),1);
        Cnhu(:,1) = MK_LOWER(row,2:end)';
        Cnhu(:,2) = MK_UPPER(row,2:end)';
        Cnhu(:,3) = MS(row,2:end)';
        Cnhu(:,4) = HZ_h141(row,2:end)';
        Cnhu(:,5) = HZ_hS(row,2:end)';
        Cnhu(:,6) = HZ_hL(row,2:end)';
        Cnhu(:,7) = WANG(row,2:end)';
        str=strcat('new_data/Cnhu_p_',num2str(p),'_n_',num2str(n_samples(i)));
        save(str, 'Cnhu')
    end
end




