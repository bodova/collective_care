%% Collect all necessary statistics for the data and save them into the folder output/DT_statistics
%% Only 6 example data files used here, one for each treatment. To get better statistics use more replicate experiments
%% Plots an ethogram for each treatment/replicate in a separate figure
%% First set the Fpath to the folder where the code is located

clear; 
global DTL DTA idty load_init typeY1
global yP yR yL NP NR NL
global Fpath

global transitions times transitionsT timesT
global loadsOWN loadsSEEN R Lreset Pactivity Ractivity
global DT ii jj VrateMM KrateMM ANTtypes prod FIX

Fpath = 'folder_path';  %% set this to the path where the folder ModelInference_Simulation is stored on the computer
run_get_stats = 1;      %% set equal to 1 to run statistics for the given memory parameters
run_infer_model = 1;    %% set equal to 1 to infer the model with the memory parameters selected in the first part
run_simulate = 1;       %% set equal to 1 to run stochastic simulation with the same rates as inferred (need to run inference at the same time)
                        %% the preset interval of simulation is 30 min, 90 min (corresponds to the experiments) runs longer and requires more patience
saveres = 1;            %% set equal to 1 to save the results into appropriate folders

%% The optimal memory parameters for the activity variable A (P,R) and for the load L (in frames = 1/15 seconds)
DTA = 15;               %% in frames = 1/15 min
DTL = 2*60*15;          %% in frames = 1/15 min


%% The mesh for all considered variables
yP = [-0.14,0.15,0.25]; NP = length(yP);
yR = [-0.14,0.15,0.25]; NR = length(yR);
yL = [0,1*10^4,5e4]; NL = length(yL);

load_init = [1,0.5,0.25,0];                  % H, L, N, M - in this order, unclassified characterized as 3rd option
repTrain = [];
repTrain{1} = [1]; repTrain{2} = [1]; repTrain{3} = [1]; repTrain{4} = [1]; repTrain{5} = [1]; repTrain{6} = [1];

%% All treatments and replicates used for the analysis
repTrain = [];
repTrain{1} = [1]; repTrain{2} = [1]; repTrain{3} = [1]; repTrain{4} = [1]; repTrain{5} = [1]; repTrain{6} = [1];
treatment = {'HH', 'HL', 'HTx', 'LL', 'LTx', 'TxTx'};

%% ----------------------------------------------------------------------------
%% ---------------> PART I: Get statistics from the experimental data
%% ----------------------------------------------------------------------------
if run_get_stats == 1
    inLfix = 1;             %% 1 when initial load set to a median experimental value, 0 initial loads set to back-calculated MM loads
    ANT = 1:6;
    TREATMENT = 1:6;
    typeY1 = 1;             % only allo activity
    
    TRANSITIONS = [];       % store transition number for all treatments
    TIMES = [];             % store occupancy time for all treatments
    fprintf(1,'Obtaining sufficient statistics:\nMemory A: %d, Memory L: %d (1s = 15, 1min = 900)',DTA,DTL)
    
    %% Universal parameters
    VrateMM = 5.22;         % Michaelis-Menten: dL = -vl/(l+K)
    KrateMM = 50119;
    type = {'POST','PRE'};
    mark = {'blue', 'red', 'green', 'purple', 'yellow', 'orange'};
    fps = 15;
    runPOST = 1;
    Lreset = 1;             % if 1 only takes last DT into consideration
    track = 0;
    dataid = num2str('1');
            
    for idtr = TREATMENT                        % index of the treatment: 'HH', 'HL', 'HTx', 'LL', 'LTx', 'TxTx'
        fprintf(1,'\nTreatment %d: ',idtr)
        for idrep = 1%:numel(repTrain{idtr})    % through all replicates
            fprintf(1,'replicate %d ',repTrain{idtr}(idrep))
            
            %% ------------> POST phase
            if runPOST == 1
                idty = 1; % POST phase
                [events,durations,treat,loadin,idallofrom,idalloto,idself] = readdish(repTrain,idty,idtr,idrep);
                [charACT_Pa,charACT_Ra,charACT_Pas,charACT_Ras,charACT_Ps,~,loadSEEN_a,~,loadSEEN_as] = get_dishstats(idallofrom,idalloto,idself,loadin,VrateMM,KrateMM);
                if typeY1 == 1
                    charACT_P = charACT_Pa;
                    charACT_R = charACT_Ra;
                elseif typeY1 == 2
                    charACT_P = charACT_Ps;
                    charACT_R = charACT_Ps;
                else
                    charACT_P = charACT_Pas;
                    charACT_R = charACT_Ras;
                end
                
                plotraster(treat,treatment{idtr},repTrain{idtr}(idrep),idallofrom,idalloto,idself);
                
                shift = 1; % 1 = no shift
                meansPOST = [];
                for ant= ANT
                    tr0 = zeros(8,1); ti0 = zeros(3,1);
                    trP = zeros(8,NP);  tiP = zeros(3,NP);
                    trR = zeros(8,NR);  tiR = zeros(3,NR);
                    trLa = zeros(8,NL); tiLa = zeros(3,NL);
                    trLas = zeros(8,NL);  tiLas = zeros(3,NL);
                    trPLa = zeros(8,NP,NL);  tiPLa = zeros(3,NP,NL);
                    trRLa = zeros(8,NR,NL); tiRLa = zeros(3,NR,NL);
                    trPLas = zeros(8,NP,NL);  tiPLas = zeros(3,NP,NL);
                    trRLas = zeros(8,NR,NL); tiRLas = zeros(3,NR,NL);
                    histPRDL = [];
                    
                    loads{ant} = [];
                    U = cellfun(@minus,charACT_P,charACT_R,'Un',0);
                    [transitions0,times0] = getstats_0D(ant,events,durations);
                    [transitionsPR,timesPR,histPR] = getstats_2D(ant,charACT_P,charACT_R,events,durations,[0,0],NP,yP,NR,yR,DTA,DTL,[0,0]);   % 2D activity P,R
                    [transitionsPLa,timesPLa,histPLa] = getstats_2D(ant,charACT_P,loadSEEN_a,events,durations,[0,0],NP,yP,NL,yL,DTA,DTL,[0,1]);   % 2D activity P,L
                    [transitionsRLa,timesRLa,histRLa] = getstats_2D(ant,charACT_R,loadSEEN_a,events,durations,[0,0],NR,yR,NL,yL,DTA,DTL,[0,1]);   % 2D activity R,L
                    [transitionsPLas,timesPLas,histPLas] = getstats_2D(ant,charACT_P,loadSEEN_as,events,durations,[0,0],NP,yP,NL,yL,DTA,DTL,[0,1]);   % 2D activity P,L
                    [transitionsRLas,timesRLas,histRLas] = getstats_2D(ant,charACT_R,loadSEEN_as,events,durations,[0,0],NR,yR,NL,yL,DTA,DTL,[0,1]);   % 2D activity R,L
                    
                    transitionsP = sum(transitionsPLas,3);
                    timesP = sum(timesPLas,3);
                    transitionsR = sum(transitionsRLas,3);
                    timesR = sum(timesRLas,3);
                    transitionsLa = reshape(sum(transitionsPLa,2),[8,3]);
                    timesLa = reshape(sum(timesPLa,2),[3,3]);
                    transitionsLas = reshape(sum(transitionsPLas,2),[8,3]);
                    timesLas = reshape(sum(timesPLas,2),[3,3]);
                    
                    tr0 = tr0 + transitions0; ti0 = ti0 + times0;
                    trP = trP + transitionsP; tiP = tiP + timesP;
                    trR = trR + transitionsR; tiR = tiR + timesR;
                    trLa = trLa + transitionsLa; tiLa = tiLa + timesLa;
                    trLas = trLas + transitionsLas; tiLas = tiLas + timesLas;
                    trPLa = trPLa + transitionsPLa; tiPLa = tiPLa + timesPLa;
                    trRLa = trRLa + transitionsRLa; tiRLa = tiRLa + timesRLa;
                    trPLas = trPLas + transitionsPLas; tiPLas = tiPLas + timesPLas;
                    trRLas = trRLas + transitionsRLas; tiRLas = tiRLas + timesRLas;
                    
                    %(treatment, replicate, ant, activity type, PRE/POST, model, load type),
                    % typeY1 = 1, 2, 3 (activity allo, self, both),
                    % loadtype = 1, 2 (load decays during allo, allo+self),
                    % model = 1=P, 2=R, 6=L, 7=PL, 8=RL, 9=C
                    % 1/2 if allo / allo+self
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,9,1} = sum(trP'); TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,9,1} = sum(tiP');
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,1,1} = trP; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,1,1} = tiP;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,2,1} = trR; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,2,1} = tiR;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,6,1} = trLa; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,6,1} = tiLa;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,6,2} = trLas; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,6,2} = tiLas;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,7,1} = trPLa; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,7,1} = tiPLa;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,8,1} = trRLa; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,8,1} = tiRLa;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,7,2} = trPLas; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,7,2} = tiPLas;
                    TRANSITIONS{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,8,2} = trRLas; TIMES{idtr,repTrain{idtr}(idrep),ant,typeY1,idty,8,2} = tiRLas;
                end
            end
        end
    end
    fprintf(1,'\n')
    if saveres == 1
        save([Fpath,'/ModelInference_Simulation/output/DT_statistics/results_DTA',num2str(DTA),'_DTL',num2str(DTL)],'TRANSITIONS','TIMES','idallofrom','idalloto','idself');
    end
end
%% ----------------------------------------------------------------------------
%% <----------------- End getting statistics (part I)
%% ----------------------------------------------------------------------------


%% ----------------------------------------------------------------------------
%% -----------------> PART II: Infer the model of our choice (need to select the model: variable variant)
%% ----------------------------------------------------------------------------

%% Model choice - need to run statistics first, uses the same memory parameters as the fisrt part of the code
if run_infer_model == 1
    fprintf(1,'\nInferring the best model')
    load([Fpath,'/ModelInference_Simulation/output/DT_statistics/results_DTA',num2str(DTA),'_DTL',num2str(DTL)],'TRANSITIONS','TIMES','idallofrom','idalloto','idself');
    
    variant = 8;        %% the best model 9=RLa, statistics already in the folder
    idty = 1;           %% POST phase
    T = 15*60*90;       %% 90 min of the experiment

    if variant == 1 %C
        mid = 5; MODin = 9; TA = 1; TL = 1; fprintf(1,'\nModel C: ')
    elseif variant == 2 %Pa
        mid = 1; MODin = 1; TA = 1; TL = 1; fprintf(1,'\nModel P: ')
    elseif variant == 3 %Ra
        mid = 2; MODin = 2; TA = 1; TL = 1; fprintf(1,'\nModel R: ')
    elseif variant == 4 %La
        mid = 4; MODin = 6; TA = 1; TL = 1; fprintf(1,'\nModel La: ')
    elseif variant == 5 %Las
        mid = 4; MODin = 6; TA = 1; TL = 2; fprintf(1,'\nModel Las: ')
    elseif variant == 6 %PaLa
        mid = [1,4]; MODin = 7; TA = 1; TL = 1; fprintf(1,'\nModel PLa: ')
    elseif variant == 7 %PaLas
        mid = [1,4]; MODin = 7; TA = 1; TL = 2; fprintf(1,'\nModel PLas: ')
    elseif variant == 8 %RaLa
        mid = [2,4]; MODin = 8; TA = 1; TL = 1; fprintf(1,'\nModel RLa: ')
    elseif variant == 9 %RaLas
        mid = [2,4]; MODin = 8; TA = 1; TL = 2; fprintf(1,'\nModel RLas: ')
    end
    
    %% ---------> KEY PARAMETERS
    TXload = 1; loadTX = 5e4;       %% load of Tx ants
    Lthresh = [5*10^4];             %% load detection threshold
    samp = 2;                       %% 1-random, 2-sequential search, 3-only highest load, nonzero detection
    %% 4-only loaded, nonzero detection, propto loadsize,
    FIX = 0;                        %% 0 = rho optimized, 1 = rho fixed
    prod = 100;                     %% rho multiplied by 100
    EPS = 0.075;                    %% parameter in the SEQ rule
    binning = 1;                    %% 1 = fine 3x3 (8)
    inLfix = 0;                     %% 0 initial load set to a back-computed MM load
    DT = max(DTA,DTL);              %% initial segment of the experimental data needed to start the simulation
    VrateMM = 5.22;                 %% decay rate of the load
    KrateMM = 50119;                %% constant in MM dynamics dL = -bv*L/(L+K) if allo groomed
    A1 = 0;                         %% maximal activity is not capped at 1
    ANTtypes1 = {'M'};              %% unmodified rates by rho
    ANTtypes2 = {'H','L','N','X'};  %% modified rates by rho
    Lreset = 1;
    
    plotload = 0;
    simulate = 0;
    infer = 1;
    plotratesE = 0;
    plotratesS = 0;
    plotActivity = 0;
    
    Atype = TA;                     %% 1=allo (Pa, Ra)
    typeL = TL;                     %% 1 = allo (La), 2 = allo + self (Las)
    rp = 1;                         %% replicate simulation
    
    MODELS = {'P','R','D','L','1'};     %% all possible 1D models
    idvar = find(mid==4,1);             %% position of load in the explanatory variable space
    if isempty(idvar)
        idvar = 1;
    end
    NVAR = length(mid);                 %% number of variables
    
    %% Set the mesh for all considered variables
    yP = [-0.14,0.15,0.25]; NP = length(yP);
    yR = [-0.14,0.15,0.25]; NR = length(yR);
    yL = [0,1*10^4,5e4]; NL = length(yL);
    
    y0 = 0;
    Ymesh = {yP,yR,yR,yL,y0};
    mesh = []; Nmesh = [];
    for i=1:NVAR
        mesh{i} = [];
        model{i} = MODELS{mid(i)};
        ifL(i) = double(model{i}=='L');
        mesh{i} = Ymesh{mid(i)};
        Nmesh(i) = length(mesh{i});
    end
    
    %% ---------> PART 2: Experimental data
    load_init = [1,0.5,0.25,0];             % H, L, N, M - in this order, unclassified characterized as 3rd option
    t = DT;                                 %initial time
    Ants = 6;                               % number of ants in one dish
    load([Fpath,'/ModelInference_Simulation/output/AntTreatments']); %% treatAnt (cell array with HLNMX types)
    
    [RATES,RHO,Lexp] = Inference_rho(idty,MODin,typeL,Atype,repTrain);
    
    trEN = transitions; tiEN = times;       %% sufficient statistics for the nestmates
    trET = transitionsT; tiET = timesT;     %% sufficient statistics for the treated ants
    rates = exp(RATES);
    if FIX == 0
        prod = RHO;
    end
    ratesEN = rates;                        %% fit on the nestmates
    ratesET = ratesEN;                      %% fit on the treated ants - initialization
    if prod > 0                             %% fit on the treated ants
        ratesET(2,:,:) = ratesEN(2,:,:)/RHO;
    end
    
    plotrates(log(ratesEN),MODin,trEN,tiEN,[0.3,0.7,0.3],7);
    %% XS, XA, SX, SA, AX, AS, SS, AA
    %% X-idle, S-self, A-allo
    fprintf(1,'\n')

end

%% ----------------------------------------------------------------------------
%% <----------------- End inferring the model (part II)
%% ----------------------------------------------------------------------------



%% ----------------------------------------------------------------------------
%% -----------------> PART III: Simulate using the inferred rates and compare with experiments
%% ----------------------------------------------------------------------------

if run_simulate == 1 
    fprintf(1,'\nRunning stochastic simulation')
    col = get(gca,'colororder');
    col = [col;[0,0,0]];
    filesid = 'RHOopt_ALL';
    variant = 9; %RaLas
    idty = 1; 
    
    T = 15*60*30;       %% length of the simulation in frames
    DTA = 1*15; DTL = 30*15; mid = [2,4]; TA = 1; TL = 2;
    
    TXload = 1;
    loadTX = 5e4;
    samp = 2;           %% 1=RAND rule, 2=SEQ rule, 3=MAX rule
    FIX = 0;            %% 0 = rho optimized, 1 = rho fixed
    prod = 100;         %% /100
    EPS = 0.075;        %% parameter in the acceptance probability (in case of a SEQ rule)
    binning = 1;        %% resolution to 3x3 rates for each transition
    inLfix = 0;         %% 0 - initial load set to the back-computed MM load
    A1 = 0;             %% 0 no constraints on the maximal activity
    
    %% ---------> KEY PARAMETERS
    DT = max(DTA,DTL);  %% need to read that much of the experimental data
    VrateMM = 5.22;     %% decay rate of the load
    KrateMM = 50119;    %% constant in MM dynamics dL = -bv*L/(L+K) if allo groomed
    ANTtypes1 = {'M'};                       %% unmodified rates
    ANTtypes2 = {'H','L','N','X'};           %% modified rates by the parameter rho
    Lreset = 1;
    
    %% changes here
    simulate = 1;
    plot_raster = 1;
    plot_load = 1;
    
    Lthresh = 5*10^4;
    Atype = TA;%[1,3]
    typeL = TL;%[1,2]
    rp = 1;                         %% replicate simulation
    idtr = 1;           %% Treatment: 1='HH', 2='HL', 3='HTx', 4='LL', 5='LTx', 6='TxTx'
    idrep = 1;          %% Replicate, range: 1:length(repTrain{idtr})
    fprintf(1,'\nTreatment: %d, Replicate: %d, ',idtr,repTrain{idtr}(idrep))
    
    MODELS = {'P','R','D','L','1'};     %% all possible 1D models
    idvar = find(mid==4,1);             %% position of load in the explanatory variable space
    if isempty(idvar)
        idvar = 1;
    end        
    NVAR = length(mid);                 %% number of variables
    
    % model number 
    if mid == 1;                MODin = 1; %Pa
    elseif mid == 2;            MODin = 2; %Ra
    elseif mid == 4;            MODin = 6; %La
    elseif mid == 5;        	MODin = 9; %0
    elseif mid == [1,4];        MODin = 7; %PaL
    elseif mid == [2,4];        MODin = 8; %RaL
    end
    
    %% Set mesh (same as for inference)
    yD = yP;
    y0 = 0;
    Ymesh = {yP,yR,yD,yL,y0};
    mesh = []; Nmesh = [];
    for i=1:NVAR
        mesh{i} = [];
        model{i} = MODELS{mid(i)};
        ifL(i) = double(model{i}=='L');
        mesh{i} = Ymesh{mid(i)};
        Nmesh(i) = length(mesh{i});
    end
    
    %% ---------> PART 2: Experimental data
    %% Read data from experiments for the initial DT
    load_init = [1,0.5,0.25,0];             % H, L, N, M - in this order, unclassified characterized as 3rd option
    t = DT;                                 %initial time
    Ants = 6;                               % number of ants in one dish
    load([Fpath,'/ModelInference_Simulation/output/AntTreatments']); %% treatAnt (cell array with HLNMX types)
    
    [events,durations,treat,loadin,idallofrom,idalloto,idself] = readdish(repTrain,idty,idtr,idrep); % Load initial DT of the experiment

    %% set loadin for Tx ants to be loadTX
    if TXload == 1
        for iant = 1:6
            if string(treat(iant)) == 'N'
                loadin{iant} = loadTX;
            end
        end
    end
    
    %% Plot the experimental ethogram
    if plot_raster == 1
        plotraster(treat,treatment{idtr},idrep,idallofrom,idalloto,idself)
        for i=1:3
            subplot(3,1,i); hold on
            line([DT,DT]/15/60,[0.5,6.5],'color','k','linewidth',1)
        end
    end
    
    idallotoE = idalloto;
    idallofromE = idallofrom;
    idselfE = idself;
    
    %% Plot the experimental loads - inferred by MM model
    loadsOWN = [];
    loadsSEEN = [];
    tim = 1:T;
    initializeL2(Atype,typeL,1,T,idallofrom,idalloto,idself,R,loadin,1); % compute initial loads
    
    loadsSEENe = loadsSEEN;
    loadsOWNe = loadsOWN;
    
    %% ---------> PART 3: Generate random ethogram following the same rules
    fprintf('\nStochastic simulation')
    % compute the initial load own/seen - to start the simulation
    loadsOWN = [];
    loadsSEEN = [];
    Pactivity = [];
    Ractivity = [];
    
    tim = 1:DT;
    S = [0,0,0,0,0,0]; % 0 = X, 1 = allo, 2 = self
    initializeL2(Atype,typeL,1,DT,idallofrom,idalloto,idself,R,loadin,1);
        
    if simulate == 1
        %% Concatenate to [0,DT]: idallofrom,idalloto,idself - initial data, further concatenate only to compute activity, not for load
        idalloto = concatenate(idalloto,0,DT);
        idallofrom = concatenate(idallofrom,0,DT);
        idself = concatenate(idself,0,DT);
        
        %% Gillespie with appending idallofrom,idalloto,idself
        S = ones(1,6);      %% 1 = X, 2 = self, 3 = allo
        for ant = 1:6
            if ~isempty(idallofrom{ant}) && idallofrom{ant}(end,3)+idallofrom{ant}(end,4) - 1 == t
                S(ant) = 3; % allo
            elseif ~isempty(idself{ant}) && idself{ant}(end,3)+idself{ant}(end,4) - 1 == t
                S(ant) = 2; % self
            else
                S(ant) = 1; % X
            end
        end
        
        % Set parameters for the stochastic simulation
        thresh = exprnd(1,1,Ants);  % exponential random variable
        K = 2*10^4; numrand = exprnd(1,K,Ants);   % generator of uniform random numbers
        iter = ones(1,Ants);                % current iteration count
        id_transition = 1;                  % number of transitions
        tim = [tim,DT+1:T];
        
        t=DT+1;
        if mid ~= 5
            Y = computeY0new(NVAR,mid,Atype,typeL,t,DT,idallofrom,idalloto,idself,loadin);
        end
        
        lOWN = cell2mat(loadin)';
        
        for t = DT+1:T %% all times
            for ant=1:Ants    % get rates and update thresh
                xAnts = 1:6; xAnts(ant) = [];
                xAntsL = 1:6;
                xAntsL(ant) = [];
                if samp<=5 %% based on updated load
                    xAntsL = find(loadsOWN(1:6,end) > Lthresh); % nonfocal ants with nonzero load
                else %% based on initial load
                    xAntsL = find(lOWN > Lthresh); % nonfocal ants with nonzero load
                end
                xAntsL(xAntsL==ant)=[];
                if ~isempty(xAntsL)
                    if samp<=5 %% based on updated load
                        xAntsLoads = loadsOWN(xAntsL,end)/sum(loadsOWN(xAntsL,end)); % weights of the infected ants proportional to their loads
                    else %% based on initial load
                        xAntsLoads = lOWN(xAntsL)/sum(lOWN(xAntsL));
                    end
                end
                if samp<=5
                    Lweights = max(loadsOWN(xAnts,end),Lthresh);
                else
                    Lweights = max(lOWN(xAnts),Lthresh);
                end
                Lweights = Lweights/sum(Lweights);
                
                
                if S(ant) == 3 % compute off rate XS, XA, SX, SA, AX, AS, SS, AA
                    Rout = [5,6,8]; % AX, AS, AA
                    Sout = [1,2,3]; % X, S, A
                elseif S(ant) == 2
                    Rout = [3,4,7]; % SX, SA, SS
                    Sout = [1,3,2]; % X, A, S
                else
                    Rout = [1,2]; % XS, XA
                    Sout = [2,3]; % S, A
                end
                
                
                
                %% use instead of rates ratesET and ratesEN for treated and nestmates need to check ant type
                if string(treatAnt{idtr,repTrain{idtr}(idrep)}(ant)) == ANTtypes1
                    rates = ratesEN; % nestmate
                else
                    rates = ratesET; % treated
                end
                
                if mid == 5
                    rates_off = rates(Rout);
                elseif NVAR == 1
                    [~, id(ant,1)] = min(abs(Y(ant,1)-mesh{1}));
                    rates_off = rates(Rout,id(ant));
                elseif NVAR == 2
                    [~, id(ant,1)] = min(abs(Y(ant,1)-mesh{1}));  % index of the Y in the mesh1
                    [~, id(ant,2)] = min(abs(Y(ant,2)-mesh{2}));  % index of the Y in the mesh2
                    rates_off = rates(Rout,id(ant,1),id(ant,2));
                end
                thresh(ant) = thresh(ant) - sum(rates_off); % rate r, time t -> exp(-r*t) = xi
                
                
                if thresh(ant) > 0 % was active, increase duration by 1
                    if S(ant) == 3 % 1=nothing, 3=allo, 2=self
                        idallofrom{ant}(end,4) = idallofrom{ant}(end,4) + 1;            % increase duration in idallofrom
                        antTO = idallofrom{ant}(end,2);
                        idalloto{antTO}(end,4) = idalloto{antTO}(end,4) + 1;                % increase duration in idalloto
                    elseif S(ant) == 2
                        idself{ant}(end,4) = idself{ant}(end,4) + 1;
                    end
                else
                    Sold = S(ant);                                                      % old state of performing ant
                    S(ant) = randsample(Sout,1,true,rates_off/sum(rates_off));          % pick a new state at random
                    % PICK RECEIVER: pick an ant to receive the new action (if not self)
                    if samp == 1
                        antTO = xAnts(randsample(5,1));                                    % Alternative 1: new receiver random
                    elseif samp == 4  || samp == 6  || samp == 8                                           % samp = 4 has nonzero load detection threshold
                        if ~isempty(xAntsL)
                            antTO = xAntsL(randsample(length(xAntsL),1,true,xAntsLoads));      % Alternative 2: new receiver nonzero load - proportional to the load size
                        else
                            antTO = xAnts(randsample(5,1));
                        end
                    elseif samp == 2    %% sequential:
                        %% xAntsL - nonfocal ants with load > Lthresh
                        %% xAntsLoads - loads of these ants
                        %% Lweights - weights, max(xAntsLoads,Lthresh)
                        if ~isempty(xAntsL)
                            %% sequential search of all ants with
                            %% acceptance probability p(L) = exp[-eps-L/Lthresh]
                            accepted = 0;
                            while accepted == 0
                                idTO = randsample(5,1);
                                antTO = xAnts(idTO); %% randomly sampled nonfocal ant
                                
                                rr = rand;
                                PA = 1-exp(-EPS-loadsOWN(antTO,end)/Lthresh);
                                if rr < PA %% same Lthresh
                                    accepted = 1;
                                end
                            end
                        else
                            antTO = xAnts(randsample(5,1));
                        end
                    elseif samp == 3 || samp == 7
                        if ~isempty(xAntsL)
                            [val,idx] = max(xAntsLoads);
                            antTO = xAntsL(idx);
                        else
                            antTO = xAnts(randsample(5,1));
                        end
                    elseif samp == 5 || samp == 9
                        if ~isempty(Lweights)
                            antTO = xAnts(randsample(5,1,true,Lweights));
                        else
                            antTO = xAnts(randsample(5,1));
                        end
                    end
                    
                    if S(ant) == 3
                        idallofrom{ant} = [idallofrom{ant}; [ant,antTO,t,1]];           % ALLO
                        idalloto{antTO} = [idalloto{antTO}; [ant,antTO,t,1]];           % ALLO
                    elseif S(ant) == 2
                        idself{ant} = [idself{ant}; [ant,ant,t,1]];      % SELF
                    end
                    thresh(ant) = numrand(iter(ant),ant);
                    iter(ant) = iter(ant)+1;
                    id_transition = id_transition + 1;
                end
            end
            
            
            [l_OWN_update,l_SEEN_update] = update_loads(S,idallofrom,loadsOWN,loadsSEEN,typeL,VrateMM,KrateMM); % current state of everyone, info on who to who
            loadsOWN = [loadsOWN,l_OWN_update'];
            loadsSEEN = [loadsSEEN,l_SEEN_update'];
            if mid ~= 5
                Y = computeYnew(loadsOWN,loadsSEEN,NVAR,mid,Atype,typeL,t,DT,idallofrom,idalloto,idself,loadin);
            end
        end
        
        if plot_raster == 1
            plotraster(treat,treatment{idtr},idrep,idallofrom,idalloto,idself)
            for i=1:3
                subplot(3,1,i); hold on
                line([DT,DT]/15/60,[0.5,6.5],'color','k','linewidth',1)
            end
            subplot(3,1,1); title('Simulation')
        end
        
        
    end
    idallotoS = idalloto;
    idallofromS = idallofrom;
    idselfS = idself;
end

%% ----------------------------------------------------------------------------
%% -----------------> End simulation (part III)
%% ----------------------------------------------------------------------------



%% ----------------------------------------------------------------------------
%% -----------------> Auxiliary functions
%% ----------------------------------------------------------------------------
function[] = plotraster(treat,TR,idrep,idallofrom,idalloto,idself)
    figure; clf; 
    ca = [100,178,100]/255; 
    
    %% performed
    subplot(3,1,1); hold on; box on; set(gca,'fontsize',14)
    title(['Treatment: ',TR,' Replicate: ',num2str(idrep)])
    ylim([0.5,6.5]); xlim([-5,90])
    yticks(1:6); set(gca,'xtick',[0,30,60,90])
    lw = 8;
    ylabel('Performed')
    
    for ant = 1:6
        IDsize = size(idallofrom{ant});
        for i=1:IDsize(1)
            line([idallofrom{ant}(i,3),idallofrom{ant}(i,3)+idallofrom{ant}(i,4)]/15/60,[ant,ant],'color',ca,'linewidth',lw)
        end
    end
    yticks([]); ms = 12; p = -2.5;
    for ant=1:6
        %ant, string(treat{ord(ant)})
        if string(treat{ant}) == 'M'  plot(p,ant,'x','markersize',ms,'color','k','linewidth',2)% nestmates
        elseif string(treat{ant}) == 'H' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor','k')
        elseif string(treat{ant}) == 'L' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[0.5,0.5,0.5])
        elseif string(treat{ant}) == 'N' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[1,1,1])
        else plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[0.5,0.5,0.5])
        end
    end

    %% received
    subplot(3,1,2); hold on; box on; set(gca,'fontsize',14)
    ylim([0.5,6.5]); xlim([-5,90])
    yticks(1:6); set(gca,'xtick',[0,30,60,90])
    ylabel('Received')
    for ant = 1:6
        IDsize = size(idalloto{ant});
        for i=1:IDsize(1)
            line([idalloto{ant}(i,3),idalloto{ant}(i,3)+idalloto{ant}(i,4)]/15/60,[ant,ant],'color',ca,'linewidth',lw)
        end
    end
    yticks([]); ms = 12; p = -2.5;
    for ant=1:6
        if string(treat{ant}) == 'M'  plot(p,ant,'x','markersize',ms,'color','k','linewidth',2)% nestmates
        elseif string(treat{ant}) == 'H' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor','k')
        elseif string(treat{ant}) == 'L' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[0.5,0.5,0.5])
        elseif string(treat{ant}) == 'N' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',0.7*[1,1,1])
        else plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[0.5,0.5,0.5])
        end
    end
  
    %% self
    subplot(3,1,3); hold on; box on; set(gca,'fontsize',14)
    ylim([0.5,6.5]); xlim([-5,90])
    yticks(1:6); set(gca,'xtick',[0,30,60,90])
    ylabel('Self')
    xlabel('Time (min)')
    for ant = 1:6
        IDsize = size(idself{ant});
        for i=1:IDsize(1)
            line([idself{ant}(i,3),idself{ant}(i,3)+idself{ant}(i,4)]/15/60,[ant,ant],'color',ca,'linewidth',lw)
        end
    end
    yticks([]); ms = 12; p = -2.5;
    for ant=1:6
        if string(treat{ant}) == 'M'  plot(p,ant,'x','markersize',ms,'color','k','linewidth',2)% nestmates
        elseif string(treat{ant}) == 'H' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor','k')
        elseif string(treat{ant}) == 'L' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[0.5,0.5,0.5])
        elseif string(treat{ant}) == 'N' plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',0.7*[1,1,1])
        else plot(p,ant,'o','markersize',ms,'color','k','linewidth',1,'markerfacecolor',[0.5,0.5,0.5])
        end
    end
end

function[loadsOWN_update,loadsSEEN_update] = update_loads(S,idallofrom,loadsOWN,loadsSEEN,typeL,VrateMM,KrateMM)

for ant=1:6
    loadsOWN_update(ant) = loadsOWN(ant,end);
    loadsSEEN_update(ant) = loadsSEEN(ant,end);
end

for ant = 1:6
    %% load always decays during both allo+self but depending on the model (typeL) we may be interested in A (1) or A+S (2)  
    if S(ant) == 2 % self, need to decay the load of ant
        antTO = ant;
        loadsOWN_update(ant) = loadsOWN_update(ant) - VrateMM*loadsOWN_update(ant)/(loadsOWN_update(ant)+KrateMM); % update own load no matter what *(typeL==1);
        if typeL == 2 %% A+S
            loadsSEEN_update(ant) = loadsOWN_update(ant);
        end
    elseif S(ant) == 3 % allo, need to decrease the load of ant TO
        antTO = idallofrom{ant}(end,2);
        loadsOWN_update(antTO) = loadsOWN_update(antTO) - VrateMM*loadsOWN_update(antTO)/(loadsOWN_update(antTO)+KrateMM);
        loadsSEEN_update(ant) = loadsOWN_update(antTO);
    end
end

end

function[] = initializeL2(Atype,Ltype,t0,t1,idallofrom,idalloto,idself,R,loadin,idseen)
global loadsOWN loadsSEEN VrateMM KrateMM Lreset Pactivity Ractivity DT DTA DTL

%% idvar is the index of the load variable - default is 1
%% Atype = [1,3], allo, allo+self
%% Ltype = [1,2], allo, allo+self

idallotoC = concatenate(idalloto,0,t1);
for ant = 1:6
    if ~isempty(idallotoC{ant})
        charallo{ant} = [];
        for i=1:length(idallotoC{ant}(:,3))
            charallo{ant} = [charallo{ant},idallotoC{ant}(:,3):idallotoC{ant}(:,3)+idallotoC{ant}(:,4)-1];
        end
    end
end

for ant=1:6
    loads_add(ant) = loadin{ant}; % define as initial load, otherwise last load - previous state
end

Pactivity = [];
Ractivity = [];
for t = t0:t1
    % compute own load and add to already computed one
    idallotoC = concatenate(idalloto,0,t);
    idallofromC = concatenate(idallofrom,0,t);
    idselfC = concatenate(idself,0,t);
    
    %% Michaelis-Menten dynamics: need to recompute at every time, based on loadin, idalloto, idallofrom, idself
    for ant = 1:6
        if ~isempty(idallotoC{ant}) && idallotoC{ant}(end,3)+idallotoC{ant}(end,4)-1 == t % at time t some ant allogrooms focal ant
            loads_add(ant) = loads_add(ant) - VrateMM*loads_add(ant)/(loads_add(ant)+KrateMM);% always update load %%*(Ltype==1||Ltype==2);
        elseif ~isempty(idselfC{ant}) && idselfC{ant}(end,3)+idselfC{ant}(end,4)-1 == t % at time t focal ant selfgrooms
            loads_add(ant) = loads_add(ant) - VrateMM*loads_add(ant)/(loads_add(ant)+KrateMM); % always update load %%*(Ltype==2); 
        end
    end
    loadsOWN = [loadsOWN,loads_add'];
    
    
    Pactivity_add = zeros(1,6);
    Ractivity_add = zeros(1,6);
    %% only the events in the last DTA
    idallotoC2 = concatenate(idalloto,max(0,t-DTA),t);
    idallofromC2 = concatenate(idallofrom,max(0,t-DTA),t);
    idselfC2 = concatenate(idself,max(0,t-DTA),t);
    dt = min(t,DTA);
    
    for ant = 1:6
        
        if ~isempty(idallofromC2{ant})
            Pactivity_add(ant) = Pactivity_add(ant) + sum(idallofromC2{ant}(:,4))/dt*(Atype==1||Atype==3);
        end
        if ~isempty(idallotoC2{ant})
            Ractivity_add(ant) = Ractivity_add(ant) + sum(idallotoC2{ant}(:,4))/dt*(Atype==1||Atype==3);
        end
        if ~isempty(idselfC2{ant})
            Pactivity_add(ant) = Pactivity_add(ant) + sum(idselfC2{ant}(:,4))/dt*(Atype==2||Atype==3);
            Ractivity_add(ant) = Ractivity_add(ant) + sum(idselfC2{ant}(:,4))/dt*(Atype==2||Atype==3);
        end
    end
    
    Pactivity = [Pactivity,Pactivity_add'];
    Ractivity = [Ractivity,Ractivity_add'];
    
    
    %% only the events in the last DTL
    idallofromC2 = concatenate(idallofrom,max(0,t-DTL),t);
    idselfC2 = concatenate(idself,max(0,t-DTL),t);
    if idseen == 1
        % compute load seen
        l_last = zeros(6,1);
        for ant = 1:6
            l_last(ant) = 0;
            %% Perceived when allo
            if Ltype == 1 %% allo
                idC = idallofromC2{ant};
            elseif Ltype == 2 %% allo+self
                if ~isempty([idallofromC2{ant};idselfC2{ant}])
                    idC = sortrows([idallofromC2{ant};idselfC2{ant}],3);
                else
                    idC = [];
                end
            end
            s = size(idC);
            if s(1)>0
                te_allo = idC(:,3) + idC(:,4);  % ending times of allo events
                a_allo = idC(:,2);                          % ants receiving allo
                SL = size(loadsOWN);
                if length(a_allo) == 1
                    l_allo = loadsOWN(a_allo,te_allo-1);
                    l_last(ant) = l_allo;
                elseif SL(2)>=te_allo-1
                    l_allo = diag(loadsOWN(a_allo,te_allo-1));                % last seen load of recent allo
                    if max(l_allo)>0
                        l_last(ant) = l_allo(find(l_allo,1,'last'));                 % value of last nonzero seen load
                    end
                end
            end
                        
        end
        loadsSEEN = [loadsSEEN,l_last]; % your own load
    end
end
end

%% Functions that compute explanatory variables at time t based on idalloto, idallofrom, idself, DTA, DTL incorporated
function[Yout] = computeY0new(NVAR,mid,Atype,Ltype,t,DT,idallofrom,idalloto,idself,loadin)
    global loadsOWN loadsSEEN Pactivity Ractivity
    global R VrateMM KrateMM Lreset
    global DTA DTL

    idallotoC = concatenate(idalloto,t-DTA,t);
    idallofromC = concatenate(idallofrom,t-DTA,t);
    idselfC = concatenate(idself,t-DTA,t);

    loads_add = loadsOWN(:,end);
    Y = zeros(6,5);
    for ant = 1:6
        if ~isempty(idallofromC{ant})
            Y(ant,1) = Y(ant,1) + ismember(Atype,[1,3])*sum(idallofromC{ant}(:,4))/DTA; % performing A
        end
        if ~isempty(idselfC{ant})
            Y(ant,1) = Y(ant,1) + ismember(Atype,[2,3])*sum(idselfC{ant}(:,4))/DTA; % performing S
        end

        if ~isempty(idallotoC{ant})
            Y(ant,2) = Y(ant,2) + ismember(Atype,[1,3])*sum(idallotoC{ant}(:,4))/DTA; % receiving A
        end
        if ~isempty(idselfC{ant})
            Y(ant,2) = Y(ant,2) + ismember(Atype,[2,3])*sum(idselfC{ant}(:,4))/DTA; % receiving S
        end

        Pactivity_add(ant) = Y(ant,1);  %P
        Ractivity_add(ant) = Y(ant,2);  %R
        Y(ant,3) = Y(ant,1) - Y(ant,2); %D

        idallotoC = concatenate(idalloto,t-DTL,t);
        idallofromC = concatenate(idallofrom,t-DTL,t);
        idselfC = concatenate(idself,t-DTL,t);


        %% Michaelis-Menten model
        if ~isempty(idallotoC{ant}) && idallotoC{ant}(end,3)+idallotoC{ant}(end,4)-1 == t % at time t some ant allogrooms focal ant
            loads_add(ant) = loads_add(ant) - VrateMM*loads_add(ant)/(loads_add(ant)+KrateMM)*(Ltype==1||Ltype==2);
        elseif ~isempty(idselfC{ant}) && idselfC{ant}(end,3)+idselfC{ant}(end,4)-1 == t % at time t focal ant selfgrooms
            loads_add(ant) = loads_add(ant) - VrateMM*loads_add(ant)/(loads_add(ant)+KrateMM)*(Ltype==2);
        end


    end
    loadsOWN = [loadsOWN,loads_add];
    Pactivity = [Pactivity,Pactivity_add'];
    Ractivity = [Ractivity,Ractivity_add'];


    % compute load seen, if Lreset = 1 only computes with last DT worth of data
    if Lreset == 0 
        IDa = idallofrom; IDs = idself;
    else
        IDa = idallofromC; IDs = idselfC;
    end

    l_last = zeros(6,1);
    t_last = 0;
    for ant = 1:6 
        if Ltype == 1 %% allo
            ID = IDa{ant};
        elseif Ltype == 2 %self
            if ~isempty([IDa{ant};IDs{ant}])
                ID = sortrows([IDa{ant};IDs{ant}],3);
            else 
                ID = [];
            end
        end
        l_last(ant) = 0;
        s = size(ID);
        if s(1)>0
            te_allo = ID(:,3) + ID(:,4);  % ending times of allo events
            a_allo = ID(:,2);                          % ants receiving allo
            SL = size(loadsOWN);
            if SL(2)>=te_allo-1
                if length(a_allo) == 1
                    l_allo = loadsOWN(a_allo,te_allo-1);
                    l_last(ant) = l_allo;
                    t_last = te_allo;
                else
                    l_allo = diag(loadsOWN(a_allo,te_allo-1));                % last seen load of recent allo
                    if max(l_allo)>0
                        l_last(ant) = l_allo(find(l_allo,1,'last'));                 % value of last nonzero seen load
                        idl = find(l_allo == l_last(ant));                           % index of the last nonzero seen load
                        idl = idl(end);
                        t_last = te_allo(idl);                                  % time of last nonzero load seen during allo
                    end
                end
            end
        end

        Y(ant,4) = l_last(ant); % last seen nonzero load*
    end
    loadsSEEN = [loadsSEEN,l_last]; % last seen nonzero load*
    for idvar=1:NVAR
        Yout(:,idvar) = Y(:,mid(idvar));
    end
end

%% Functions that compute explanatory variables at time t based on idalloto, idallofrom, idself, DTA, DTL incorporated
function[Yout] = computeYnew(loadsOWN,loadsSEEN,NVAR,mid,Atype,Ltype,t,DT,idallofrom,idalloto,idself,loadin)
global  Lreset DTA DTL
idallotoC = concatenate(idalloto,t-DTA,t);
idallofromC = concatenate(idallofrom,t-DTA,t);
idselfC = concatenate(idself,t-DTA,t);

Y = zeros(6,5);
for ant = 1:6
    if ~isempty(idallofromC{ant})
        Y(ant,1) = Y(ant,1) + ismember(Atype,[1,3])*sum(idallofromC{ant}(:,4))/DTA; % performing A
    end
    if ~isempty(idselfC{ant})
        Y(ant,1) = Y(ant,1) + ismember(Atype,[2,3])*sum(idselfC{ant}(:,4))/DTA; % performing S
    end
    
    if ~isempty(idallotoC{ant})
        Y(ant,2) = Y(ant,2) + ismember(Atype,[1,3])*sum(idallotoC{ant}(:,4))/DTA; % receiving A
    end
    if ~isempty(idselfC{ant})
        Y(ant,2) = Y(ant,2) + ismember(Atype,[2,3])*sum(idselfC{ant}(:,4))/DTA; % receiving S
    end
    Y(ant,3) = Y(ant,1) - Y(ant,2); %D
end


% compute last nonzero load seen, if Lreset = 1 only computes with last DT worth of data
if Lreset == 0 ID = loadsSEEN;
else ID = loadsSEEN(:,end-DTL:end);
end

for ant = 1:6
    id = ID(ant,:);
    if max(id) == 0
        l_last = 0;
    else
        l_last = id(find(id,1,'last'));
    end
    Y(ant,4) = l_last; % last seen nonzero load*
end

for idvar=1:NVAR
    Yout(:,idvar) = Y(:,mid(idvar));
end
%Yout(:,idvar) = Y(:,mid(idvar));
end

%% Function that concatenates the experimental data to interval [0,DT]
function[IDout] = concatenate(IDin,DT0,DT)
ID = IDin;
for ant=1:6
    if ~isempty(ID{ant})
        ID{ant} = ID{ant}(ID{ant}(:,3)+ID{ant}(:,4) >= DT0,:);
        ID{ant} = ID{ant}(ID{ant}(:,3) <= DT,:);
        
        idS = find(ID{ant}(:,3)<DT0);
        idE = find(ID{ant}(:,3)+ID{ant}(:,4)>DT);
        
        for i=1:length(idS)
            if ID{ant}(idS(i),3) < DT0
                ID{ant}(idS(i),4) = ID{ant}(idS(i),3) + ID{ant}(idS(i),4) - DT0 + 1;
                ID{ant}(idS(i),3) = DT0;
            end
        end
        for i=1:length(idE)
            if ID{ant}(idE(i),3) + ID{ant}(idE(i),4) > DT
                ID{ant}(idE(i),4) = DT - ID{ant}(idE(i),3) + 1;
            end
        end
    end
end
IDout = ID;
end



