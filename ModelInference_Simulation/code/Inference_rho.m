function[RATES,RHO,LTRAIN] = Inference_model_rho_yoga(TYin,MODin,ldallo,typeY1,repTrain);
%% WORK IN PROGRESS
%% Model: states X-idle, S-self, A-allo
%% Transitions: XS, XA, SX, SA, AX, AS, SS, AA

global DTA DTL idty ii jj N N1 N2
global yP yR NR NP NL ND
global times transitions timesT transitionsT
global penalization opts ANTtypes
global prod FIX Fpath

Model = MODin;
idty = TYin;

%% OPTIONS
penalization = 0;
fglik = 2;
traintest = 0;                          % 0 - use all data, 1 - use train/test splits
replicates = 1;                         % will be stored later in Nsplits: number of replicates in the train/test splitting
plotinf = 0;                            % plot transitions # and rates
plotrates = 1;                          % plot all rates for EVERY trainset
plotlikdt = 1;
randcol = 0.8*[1,1,1];
set_shift = 0;
opts = 1; %1=add, 0=mult
if FIX == 1
    Rin = log(prod/100);
else
    Rin = 0;
end

LTRAIN = []; LTEST = []; RATES = [];
LTRAIN0 = []; LTEST0 = []; RATES0 = [];
load([Fpath,'/ModelInference_Simulation/output/DT_statistics/results_DTA',num2str(DTA),'_DTL',num2str(DTL)]); % -2,-1,0,1,2,3


%% SET PARAMETERS:
type = {'POST','PRE'};                  % phase
typeY2 = ldallo;                        % 1 = allo activity, 2 = self activity, 3 = allo+self activity
loadallo = typeY2;                      % 1 = load decays during allo, 2 = allo+self load decay (this refers to when load decays, not when it is perceived. Perceived only when allo)
TR = [1:6]; Ntr = length(TR);           % treatments that go into inference
repfrom = 1; repto = 100;               % replicates to include - all (repfrom = 1; repto = 100;)
ANTtypes = {'H','L','N','M','X'};       % ant types that go into inference (H-high, L-low, N-no, M-nestmates)
ANTtypes1 = {'M'};                      % inferred rates
ANTtypes2 = {'H','L','N','X'};          % modified rates

if idty == 2 iter = [9,1:4,20:22];      % models for the PRE phase
else iter = [9,1:8,20:22];              % models for the POST phase [1:9,11:15,17:18,20:22];
end
iterL = [5:8,15,17:18];                 % models that include perturbation load L
load([Fpath,'/ModelInference_Simulation/output/AntTreatments']);  % which ant has which treatment - used in selecting ants of interest

%% model
typeModel = Model;
options = load_libraries;
sh = 0;


%% COMPUTATION FOR ALL DATA
RATES = [];
NP = length(yR); NR = NP;
TransitionsP = zeros(8,NP);     TimesP = zeros(3,NP);
TransitionsR = zeros(8,NR);     TimesR = zeros(3,NR);
TransitionsL = zeros(8,NL);     TimesL = zeros(3,NL);
TransitionsPL = zeros(8,NP,NL); TimesPL = zeros(3,NP,NL);
TransitionsRL = zeros(8,NR,NL); TimesRL = zeros(3,NR,NL);

Npositive1 = 0;  % number of ants that are of the type ANTtypes1
Npositive2 = 0;  % number of ants that are of the type ANTtypes2
Nall = 0;        % number of all ants

if double(ismember(Model,[5:8,15,17:18])) == 1 ld = loadallo; % load decays during 1=A, 2=A+S
else ld = 1;
end

transitions1 = 0*TRANSITIONS{1,1,1,1,idty,typeModel,ld}; times1 = 0*TIMES{1,1,1,1,idty,typeModel,ld};
transitions2 = 0*TRANSITIONS{1,1,1,1,idty,typeModel,ld}; times2 = 0*TIMES{1,1,1,1,idty,typeModel,ld};

%% Read data for the current split - only for ANTtypes (now only for nestmates)
for idtr = TR   %% ITERATION TREATMENTS
    typeAS = typeY1;
    for idtypeY1 = typeAS;
        for idrep = 1:length(repTrain{idtr}) %% ITERATION REPLICATES
            inSET1 = double(ismember(treatAnt{idtr,repTrain{idtr}(idrep)},ANTtypes1));      % checks if the current ant has a correct type in ANTtypes
            inSET2 = double(ismember(treatAnt{idtr,repTrain{idtr}(idrep)},ANTtypes2));      % checks if the current ant has a correct type in ANTtypes
            Npositive1 = Npositive1 + sum(inSET1);                                          % traces number of included ants
            Npositive2 = Npositive2 + sum(inSET2);
            for ant = find(inSET1 == 1)
                transitions1 = transitions1 + TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,idtypeY1,idty,typeModel,ld};
                times1 = times1 + TIMES{idtr,repTrain{idtr}(idrep),ant,idtypeY1,idty,typeModel,ld};
            end
            for ant = find(inSET2 == 1)
                transitions2 = transitions2 + TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,idtypeY1,idty,typeModel,ld};
                times2 = times2 + TIMES{idtr,repTrain{idtr}(idrep),ant,idtypeY1,idty,typeModel,ld};
            end
        end
    end
end

% Choice of the model
if (typeModel == 1) || (typeModel == 11) N = NP;
elseif (typeModel == 2) || (typeModel == 12) N = NR;
elseif (typeModel == 3) || (typeModel == 13) N1 = NP; N2 = NR;
elseif (typeModel == 4) || (typeModel == 14) N = ND;
elseif (typeModel == 5) || (typeModel == 15) N1 = ND; N2 = NL;
elseif (typeModel == 7) || (typeModel == 17) N1 = NP; N2 = NL;
elseif (typeModel == 8) || (typeModel == 18) N1 = NR; N2 = NL;
elseif typeModel == 6 N = NL;
elseif typeModel == 9
elseif (typeModel == 20) N = NP;
elseif (typeModel == 21) N1 = NP; N2 = NP;
elseif (typeModel == 22) N1 = NR; N2 = NP;
end


%% TRAINING OF THE MODEL
defval = -20;
Ltrain = 0;
collectL = [];

if typeModel == 9 % constant model does not have R
    transitions = transitions1+transitions2; times = times1+times2;
    a_separate1 = zeros(8,1);
    c0 = ones(1,8)';
    c = minFunc(@lik0D, c0,options);
    for kk=1:8
        if transitions(kk) == 0
            c(kk) = defval;
        end
    end
    Ltrain1 = -lik0D(c);
    Ltrain = Ltrain1/(Npositive1+Npositive2);
    
    collectL = [collectL; Ltrain];
    a_separate1 = c;
    a_separate2 = c;
    if set_shift == 1
        Lshift = Ltrain;
    end
    Rinf = 0;
elseif (typeModel == 1) || (typeModel == 2) || (typeModel == 4) || (typeModel == 6) 
    a_separate1 = zeros(8,N);
    a_separate2 = zeros(8,N);
    transitions = transitions1; times = times1;
    transitionsT = transitions2; timesT = times2;
    
    c0 = ones(1,8*N);
    c0(8*N+1) = 0; %% initial rho additive
    c0(end)=Rin;
    clin = minFunc(@lik1D,c0,options); 
    
    % conversion
    c = zeros(8,N);
    for i=1:8
        c(i,:) = clin(1+(i-1)*N:i*N);
    end
    Rinf = clin(end);
    
    %c = reshape(clin,N,8+1)';
    for kk=1:8
        for ii=1:numel(transitions(kk,:))
            if transitions(kk,ii)+transitionsT(kk,ii) == 0
                c(kk,ii) = defval;
            end
        end
    end
    a_separate1 = c;
    a_separate2 = a_separate1;
    
    %opt = 0; %1=add, 0=mult
    if opts == 0
        a_separate2(2,:) = a_separate1(2,:)/Rinf;
    else
        a_separate2(2,:) = a_separate1(2,:) - Rinf;
    end
    
    %conversion
    clin = [];
    for i=1:8
        clin = [clin,c(i,:)];
    end
    clin = [clin,Rinf];
    
    Linf = -lik1D(clin);
    Ltrain = Linf/(Npositive1+Npositive2);
    collectL = [collectL; Linf];
else
    % separate inference 2D
    transitions = transitions1+transitions2;
    times = times1+times2;
    a_separate = zeros(8,N1,N2);
    %fprintf(1,'\n')
    Ltmp = 0;
    for ii = 1:N1
        for jj = 1:N2
            c0 = ones(1,8)';
            c = minFunc(@lik2D, c0,options);
            Ltmp = Ltmp - lik2D(c);
            a_separate(:,ii,jj) = c;
        end
    end
    
    a_separate1 = zeros(8,N1,N2);
    a_separate2 = zeros(8,N1,N2);
    transitions = transitions1; times = times1;
    transitionsT = transitions2; timesT = times2;
    alin = shapeMtoV(a_separate,8,N1,N2);
    Ltmp = -lik2DR([alin,Rin])/(Npositive1+Npositive2);
    
    
    c0 = -5*ones(8,N1,N2);
    c0 = shapeMtoV(c0,8,N1,N2)';
    c0 = [c0;Rin];
    clin = minFunc(@lik2DR,c0,options); % all X->A transitions altered
    
    Rinf = clin(end);
    cmat = shapeVtoM(clin,8,N1,N2);
    for kk=1:8
        for ii=1:N1
            for jj=1:N2
                if transitions(kk,ii,jj)+transitionsT(kk,ii,jj) == 0
                    cmat(kk,ii,jj) = defval;
                end
            end
        end
    end
    
    a_separate1 = cmat;
    a_separate2 = a_separate1;
    
    if opts == 0
        a_separate2(2,:,:) = a_separate1(2,:,:)/Rinf;
    else
        a_separate2(2,:,:) = a_separate1(2,:,:) - Rinf;
    end
    clin = shapeMtoV(cmat,8,N1,N2)';
    clin = [clin;Rinf];
    Linf = -lik2DR(clin);
    Ltrain = Linf/(Npositive1+Npositive2);
    collectL = [collectL; Linf];
end
fprintf(1,'Ltrain = %1.8e, R = %1.2f, exp(R) = %1.2f',Ltrain,Rinf,exp(Rinf))
ratesA1 = a_separate1;
ratesA2 = a_separate2;
LTRAIN = [LTRAIN,Ltrain];
RATES = ratesA1;
RHO = exp(Rinf);



