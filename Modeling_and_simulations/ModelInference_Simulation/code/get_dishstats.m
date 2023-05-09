function[charACT_Pa,charACT_Ra,charACT_Pas,charACT_Ras,charACT_Ps,load_as,loadSEEN_a,loadSEEN_s,loadSEEN_as] = ...
    get_dishstats(idallofrom,idalloto,idself,loadin,VrateMM,KrateMM);

treatment = {'HH', 'HL', 'HTx', 'LL', 'LTx', 'TxTx'};
mark = {'blue', 'red', 'green', 'purple', 'yellow', 'orange'};
fps = 15;
track = 0;
charACT_Pa = []; charACT_Ps = [];       %% providing activity 
charACT_Ra = [];                        %% receiving activity
charACT_Pas = []; charACT_Ras = []; 
track_in = track;

%% Maximum number of frames in the current replicate
MaxSize = 0;
for idmark = 1:6
    j = track_in+idmark;        % id of the current ant
    if length(idalloto{j})>0
        MaxSize = max(max(idalloto{j}(:,3)+idalloto{j}(:,4)),MaxSize); % last active frame in the dish
    end
end
for idmark = 1:6
    j = track_in+idmark;        % id of the current ant
    if length(idself{j})>0
        MaxSize = max(max(idself{j}(:,3)+idself{j}(:,4)),MaxSize); % last active frame in the dish
    end
end

%% Computes activity as a 1/0 array
for idmark = 1:6
    j = track_in+idmark;
    charACT_Pa{j} = zeros(MaxSize,1);           % activity allo performed
    charACT_Ra{j} = zeros(MaxSize,1);           % activity allo received
    charACT_Ps{j} = zeros(MaxSize,1);           % activity self performed (= received)
    charACT_Pas{j} = zeros(MaxSize,1);          % activity allo+self performed
    charACT_Ras{j} = zeros(MaxSize,1);          % activity allo+self received
    
    IDs = idallofrom{j}; [rows,cols] = size(IDs);      % rows - number of active allo events provided by the current ant
    for irow = 1:rows
        tin = IDs(irow,3)+1;
        tout = IDs(irow,3)+IDs(irow,4);
        charACT_Pa{j}(tin:tout) = charACT_Pa{j}(tin:tout) + 1;
        charACT_Pas{j}(tin:tout) = charACT_Pas{j}(tin:tout) + 1;
    end
    
    IDs = idalloto{j}; [rows,cols] = size(IDs);        % rows - number of received allo events received by the current ant
    for irow = 1:rows
        tin = IDs(irow,3)+1;
        tout = IDs(irow,3)+IDs(irow,4);
        charACT_Ra{j}(tin:tout) = charACT_Ra{j}(tin:tout) + 1;
        charACT_Ras{j}(tin:tout) = charACT_Ras{j}(tin:tout) + 1;
    end
    
    IDs = idself{j}; [rows,cols] = size(IDs);        % rows - number of received allo events received by the current ant
    for irow = 1:rows
        tin = IDs(irow,3)+1;
        tout = IDs(irow,3)+IDs(irow,4);
        charACT_Ps{j}(tin:tout) = charACT_Ps{j}(tin:tout) + 1;
        charACT_Pas{j}(tin:tout) = charACT_Pas{j}(tin:tout) + 1;
        charACT_Ras{j}(tin:tout) = charACT_Ras{j}(tin:tout) + 1;
    end
end


%% computes perturbation load
load_as = [];       % own load - decays with allo & self 
loadSEEN_a = [];    % seen load - only seen during allo
loadSEEN_s = [];    % seen load - only seen during self
loadSEEN_as = [];   % seen load - seen during allo & self 

%% Computes perturbation load of each ant based on interaction
track_in = 0;
for idmark = 1:6
    j = track_in+idmark;
    load_as{j} = zeros(MaxSize,1);
    
    % Michaelis-Menten
    load_as{j}(1) = loadin{j};
    
    for it = 2:MaxSize % need to iterate the MM dynamics, time units are frames
        load_as{j}(it) = load_as{j}(it-1) - VrateMM*load_as{j}(it-1)/(load_as{j}(it-1)+KrateMM)*charACT_Ras{j}(it);
    end
end

%% Computes experienced perturbation load - what each ant sees
for idmark = 1:6
    j = track_in+idmark;
    loadSEEN_a{j} = zeros(MaxSize,1);
    loadSEEN_as{j} = zeros(MaxSize,1);
    loadSEEN_s{j} = zeros(MaxSize,1);
    
    IDs = idallofrom{j}; [rows,cols] = size(IDs);        % rows - number of received allo events received by the current ant
    for irow = 1:rows
        tin = IDs(irow,3)+1;
        tout = IDs(irow,3)+IDs(irow,4);
        loadSEEN_a{j}(tin:tout) = load_as{ IDs(irow,2) }(tin:tout);
        loadSEEN_as{j}(tin:tout) = load_as{ IDs(irow,2) }(tin:tout);
    end
    
    IDs = idself{j}; [rows,cols] = size(IDs);        % rows - number of received self events received by the current ant
    for irow = 1:rows
        tin = IDs(irow,3)+1;
        tout = IDs(irow,3)+IDs(irow,4);
        loadSEEN_as{j}(tin:tout) = load_as{j}(tin:tout);
        loadSEEN_s{j}(tin:tout) = load_as{j}(tin:tout);
    end
end


