function [Lik,LikGrad] = lik0D(a);
% inference of individual parameters (one-by-one) 
% determined by indices ii
% transitions rates: exp(a)

global times transitions
idx = [1,1,2,2,3,3,2,3];
Lik = 0;

for k=1:8
    Lik = Lik + a(k)*sum(transitions(k)) - exp(a(k))*sum(times(idx(k)));
end

LikGrad = zeros(1,8);
for k=1:8
    LikGrad(k) = sum(transitions(k)) - exp(a(k))*sum(times(idx(k)));
end

Lik = -Lik;
LikGrad = -LikGrad';




