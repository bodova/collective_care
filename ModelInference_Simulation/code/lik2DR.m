function [Lik,LikG] = lik2DR(alin);
% inference of individual parameters (all together because R=a(end)=a(9) is a global parameter) 
% transitions rates: exp(a)
% times, transitions - data for nestmates
% timesT, transitionsT - data for treated ants
% alin has dimension 1xN1*N2*8+1 and rates of each rtansition types form
% the first N1xN2, second, ... entries. The last entry is R

global times transitions timesT transitionsT
global N1 N2 opts
global penalization FIX

idx = [1,1,2,2,3,3,2,3];
Lik = 0;
LikGrad = zeros(8,N1,N2);
LikGradR = 0;
delta = penalization; %penalization for nonzero entries

R = alin(end);
a = shapeVtoM(alin,8,N1,N2);
%opt = 1; %1=add, 0=mult

for i = 1:N1
    for j = 1:N2
        Ltmp = 0;
        for k=1:8
            %% nestmates
            Lik = Lik + a(k,i,j)*transitions(k,i,j) - exp(a(k,i,j))*times(idx(k),i,j);
            LikGrad(k,i,j) = transitions(k,i,j) - exp(a(k,i,j))*times(idx(k),i,j);
            Ltmp = Ltmp +  a(k,i,j)*transitions(k,i,j) - exp(a(k,i,j))*times(idx(k),i,j);
            
            Lik = Lik - delta*exp(a(k,i,j));
            LikGrad(k,i,j) = LikGrad(k,i,j) - delta*exp(a(k,i,j));
        end
        
        %% treated ants
        for k=[1,3:8]
            Lik = Lik + a(k,i,j)*transitionsT(k,i,j) - exp(a(k,i,j))*timesT(idx(k),i,j);
            Ltmp = Ltmp + a(k,i,j)*transitionsT(k,i,j) - exp(a(k,i,j))*timesT(idx(k),i,j);
            LikGrad(k,i,j) = LikGrad(k,i,j) + transitionsT(k,i,j) - exp(a(k,i,j))*timesT(idx(k),i,j);
        end
        
        if opts == 0
            Lik = Lik + a(2,i,j)/alin(end)*transitionsT(2,i,j) - exp(a(2,i,j)/alin(end))*timesT(idx(2),i,j);
            LikGrad(2,i,j) = LikGrad(2,i,j) + transitionsT(2,i,j)/alin(end) - exp(a(2,i,j)/alin(end))/alin(end)*timesT(idx(2),i,j);
        end
        if opts == 1
            Lik = Lik + (a(2,i,j)-alin(end))*transitionsT(2,i,j) - exp(a(2,i,j)-alin(end))*timesT(idx(2),i,j);
            LikGrad(2,i,j) = LikGrad(2,i,j) + transitionsT(2,i,j) - exp(a(2,i,j)-alin(end))*timesT(idx(2),i,j);
        end
            
            
        
        if opts == 0
            LikGradR = LikGradR - a(2,i,j)/(alin(end)^2)*transitionsT(2,i,j) + a(2,i,j)/(alin(end)^2)*exp(a(2,i,j)/alin(end))*timesT(idx(2),i,j);
        else
            LikGradR = LikGradR - transitionsT(2,i,j) + exp( a(2,i,j)-alin(end) )*timesT(idx(2),i,j);
        end
    end
end


if FIX == 1
    LikGradR = 0;
end


LikG = shapeMtoV(LikGrad,8,N1,N2);

LikG = [LikG,LikGradR];       
Lik = -Lik;
LikG = -LikG'; 


