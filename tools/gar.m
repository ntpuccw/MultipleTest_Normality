s = RandStream('mt19937ar','Seed', sum(clock));
% s = RandStream('mt19937ar');
if verLessThan('matlab','8.0')
    RandStream.setDefaultStream(s);
else
    RandStream.setGlobalStream(s);
end
mu = [2 3 4 5 3 2];
sigma = [5 1.5 0 2 0 0;
1.5 3 1 0 0 0;
0 1 5 0 0 0;
2 0 0 2 0 0;
0 0 0 0 5 0;
0 0 0 0 0 3
];
x = mvnrnd(mu,sigma,60);
sd = std(x);
x(10,1) = x(10,1) + 3*sd(1);
x(15,1) = x(15,1) - 3*sd(1);
xnew = x([1:9 11:14 16:size(x,1)],:);
xbig = mvnrnd(mu,sigma,1000);

for i=20:30
[Eta, pval] = multi_norm_SR(x,i);
[i Eta]
end