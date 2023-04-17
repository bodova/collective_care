function [Lik,LikG] = Inference_likPREactivityR_yoga(alin);
% inference of individual parameters (all together, R=a(9)=a(end) is a global parameter, need to do it together) 
% transitions rates: exp(a)
% times, transitions - data for nestmates
% timesT, transitionsT - data for treated ants

global times transitions timesT transitionsT
global N opts
global penalization FIX

idx = [1,1,2,2,3,3,2,3];
Lik = 0;
LikGrad = zeros(8,N);
LikG = [];
LikGradR = 0;
a = zeros(8,N);
for i=1:8
    a(i,:) = alin(1+(i-1)*N:i*N);
end
R = alin(end);
delta = penalization; %penalization for nonzero entries
%opt = 1; %1=add, 0=mult


for i = 1:N
    Ltmp = 0;
    %% nestmates
    for k=1:8
        Lik = Lik + a(k,i)*transitions(k,i) - exp(a(k,i))*times(idx(k),i);
        Ltmp = Ltmp + a(k,i)*transitions(k,i) - exp(a(k,i))*times(idx(k),i);
        LikGrad(k,i) = transitions(k,i) - exp(a(k,i))*times(idx(k),i);

    %% treated ants
        if k==2
            if opts == 0
            	Lik = Lik + a(k,i)/R*transitionsT(k,i) - exp(a(k,i)/R)*timesT(idx(k),i);
                Ltmp = Ltmp + a(k,i)/R*transitionsT(k,i) - exp(a(k,i)/R)*timesT(idx(k),i);
            else
                Lik = Lik + (a(k,i)-R)*transitionsT(k,i) - exp(a(k,i)-R)*timesT(idx(k),i);
                Ltmp = Ltmp + (a(k,i)-R)*transitionsT(k,i) - exp(a(k,i)-R)*timesT(idx(k),i);
            end
        else
            Lik = Lik + a(k,i)*transitionsT(k,i) - exp(a(k,i))*timesT(idx(k),i);
            Ltmp = Ltmp + a(k,i)*transitionsT(k,i) - exp(a(k,i))*timesT(idx(k),i);
        end
        if k==2
            if opts == 0
            	LikGrad(k,i) = LikGrad(k,i) + transitionsT(k,i)/R - exp(a(k,i)/R)/R*timesT(idx(k),i);
            else
                LikGrad(k,i) = LikGrad(k,i) + transitionsT(k,i) - exp(a(k,i)-R)*timesT(idx(k),i);
            end
        else
            LikGrad(k,i) = LikGrad(k,i) + transitionsT(k,i) - exp(a(k,i))*timesT(idx(k),i);
        end
    end
    
    
    if opts == 0
        LikGradR = LikGradR - a(2,i)/(R^2)*transitionsT(2,i) + a(2,i)/(R^2)*exp(a(2,i)/R)*timesT(idx(2),i);
    else
        LikGradR = LikGradR - transitionsT(2,i) + exp(a(2,i) - R)*timesT(idx(2),i);
    end
    
    for k=1:8
        Lik = Lik - delta*exp(a(k,i));
    end
    for k=1:8
        LikGrad(k,i) = LikGrad(k,i) - delta*exp(a(k,i));
    end
end

if FIX == 1
    LikGradR=0;
end

for i=1:8
    LikG = [LikG,LikGrad(i,:)];
end
LikG = [LikG,LikGradR];

Lik = -Lik;
LikG = -LikG';



