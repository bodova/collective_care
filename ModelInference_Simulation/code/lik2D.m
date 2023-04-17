function [Lik,LikGrad] = Inference_likPREactivity2D(a);
% inference of individual parameters (one-by-one) 
% determined by indices ii
% transitions rates: exp(a)

global times transitions
global ii jj
global penalization

i=ii; j=jj;
defval = -100;
idx = [1,1,2,2,3,3,2,3];
delta = penalization;
Lik = 0;
LikGrad = zeros(1,8);

for k=1:8
    Lik = Lik + a(k)*transitions(k,i,j) - exp(a(k))*times(idx(k),i,j);
    LikGrad(k) = transitions(k,i,j) - exp(a(k))*times(idx(k),i,j);
    %penalization for nonzero entries
    Lik = Lik - delta*exp(a(k));
    LikGrad(k) = LikGrad(k) - delta*exp(a(k));
end

Lik = -Lik;
LikGrad = -LikGrad';




