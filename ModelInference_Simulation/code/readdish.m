function[events,durations,treat,loadin,idallofrom,idalloto,idself] = readdish(replicates,idty,idtr,idrep)
% Reads the individual ant data for a single dish
global Fpath

type = {'POST','PRE'};
treatment = {'HH', 'HL', 'HTx', 'LL', 'LTx', 'TxTx'};
mark = {'blue', 'red', 'green', 'purple', 'yellow', 'orange'};
fps = 15;
loadallo = []; 
loadself = []; 
idallofrom = []; 
idalloto = []; 
idself = []; 
events = [];
durations = [];
treat = []; 
loadin = []; 
DT0 = 0;
track = 0;

path_full = [Fpath,'/ModelInference_Simulation/data/EVENTS_'];
filename = strjoin([path_full,type(idty),'_',treatment(idtr),'_',num2str(replicates{idtr}(idrep)),'.csv'],'');
T = importfile(filename, 1, inf);
rows = height(T);
TM = max(T.frame(:))/fps/60;    % total time of the experiment in minutes
                                % idmark = 1; % 1:blue, 2:red, 3:green, 4:purple, 5:yellow, 6:orange
track_in = track;               % starting index of this replicate




%% PERFORMED EVENTS AND DURATIONS
for idmark = 1:6    % Reads all information about SELF (only through performed activity)
    TmarkAll = sortrows(T(string(T.A) == string(mark(idmark)),:),3);    % all events of given ant (only acting!!!)
    Tstroke = TmarkAll(string(TmarkAll.action) == 'A',:);               % all antenal stroke events
    Tmark = TmarkAll(string(TmarkAll.action) ~= 'A',:);                 % all other events
    Tmark = Tmark(string(Tmark.action) ~= 'T',:);                       % all tropholaxis events excluded
    Tmark = Tmark(string(Tmark.action) ~= 'E',:);                       % all E events excluded
    Tmark = Tmark(string(Tmark.action) ~= 'R',:);                        % all gaster grooming events excluded
    
    if DT0 > 0 % read only from [0,DT0]
        Tmark = Tmark(Tmark.frame<DT0,:);
    end
    track = track + 1;
    treat{track} = [];
    events{track} = [];
    durations{track} = [];
    strokes{track} = [];    % events of antenal strokes (times)
    loadin{track} = [];     % initial condition for the perturbation load
    strokes{track} = Tstroke.frame; % antenal stroke times (duration 1 frame)
    
    %% TREAT: status of the currently investigated ant, X - unidentified
    if height(Tmark)>=1    %% check!!!
        if string(Tmark.Atreat(1)) == 'N'       % non-treated (G = green, R = red, N = non-treated)
            treat{track} = 'M';                 % nestmate
        elseif string(Tmark.Alevel(1)) == 'H'   % treated high dose
            treat{track} = 'H';
        elseif string(Tmark.Alevel(1)) == 'L'   % treated low dose
            treat{track} = 'L';
        else treat{track} = 'N';                % treated no dose
        end
    else treat{track} = 'X';
    end
    
    %% EVENTS, DURATIONS
    ST = size(Tmark);
    if ST(1) > 0
        % initial X state
        if Tmark.frame(1) > 0 % if the first action starts at time > 0 (idle before)
            events{track} = [events{track};'X']; % idle initial state
            durations{track} = [durations{track};[0,Tmark.frame(1)-1,Tmark.frame(1)]]; % until initial time
        end
        for i=1:height(Tmark)-1
            events{track} = [events{track};string(getstate(Tmark.action(i)))];
            durations{track} = [durations{track};[Tmark.frame(i),Tmark.frame(i)+Tmark.duration(i)-1,Tmark.duration(i)]];
            if Tmark.frame(i+1) > Tmark.frame(i)+Tmark.duration(i)
                events{track} = [events{track};'X'];
                durations{track} = [durations{track};[Tmark.frame(i)+Tmark.duration(i),Tmark.frame(i+1) - 1,Tmark.frame(i+1)-Tmark.frame(i)-Tmark.duration(i)]];
            end            
        end
        for i=1:size(durations{track},1)-1
            if durations{track}(i+1,1) == durations{track}(i,2)
                durations{track}(i,2) = durations{track}(i,2)-1;
                durations{track}(i,3) = durations{track}(i,3)-1;
            end
        end
        if DT0 > 0  % add the last event
            if Tmark.frame(ST(1))+Tmark.duration(ST(1)) > DT0 % event on at DT0
                durations{track} = [durations{track};[Tmark.frame(ST(1)),DT0,DT0-Tmark.frame(ST(1))+1]];
                events{track} = [events{track};string(getstate(Tmark.action(ST(1))))];
            else
                durations{track} = [durations{track};[Tmark.frame(ST(1)),Tmark.frame(ST(1))+Tmark.duration(ST(1)),Tmark.duration(ST(1))]];
                events{track} = [events{track};string(getstate(Tmark.action(ST(1))))];
                durations{track} = [durations{track};[Tmark.frame(ST(1))+Tmark.duration(ST(1))+1,DT0,DT0-Tmark.frame(ST(1))-Tmark.duration(ST(1))-1]];
                events{track} = [events{track};'X'];
            end
        end
    end
end






%% Goal: obtain for each ant the encountered perturbation load
LDin = readtable('/Users/kbodova/work/Work_Ants_Barbara/matlab2/data/initial_load_MM.csv');        % read experimental loads from the load file
LDin = LDin(string(LDin.treat) == treatment(idtr),:);               % take only the ones relevant for the current treatment
LDin = LDin(LDin.rep == replicates{idtr}(idrep),:);
[dim1,dim2] = size(LDin);

for idmark = 1:6
        loadallo{track_in+idmark} = [];     % no ALLO-events initially
        loadself{track_in+idmark} = [];     % no SELF-events initially
        idallofrom{track_in+idmark} = [];   % index of the ant making action
        idalloto{track_in+idmark} = [];     % index of the ant receiving action
        idself{track_in+idmark} = [];       % index of the self-cleaning ant
        loadin{track_in+idmark} = [0];      % initial load
    if dim1 > 0
        if idty == 1
            %string(mark{idmark})
            if string(mark{idmark}) == LDin.color{1}
                loadin{idmark} = LDin.spores(1);
            elseif string(mark{idmark}) == LDin.color{2}
                loadin{idmark} = LDin.spores(2);
            end
        else
            loadin{idmark} = 0;
        end
    end
end

for idmark = 1:6
    idfrom = idmark;
    TmarkAll = sortrows(T(string(T.A) == string(mark(idmark)),:),3);    % all events of given ant
    Tmark = TmarkAll(string(TmarkAll.action) ~= 'A',:);                 % all other events
    Tmark = Tmark(string(Tmark.action) ~= 'T',:);                       % all tropholaxis events excluded
    Tmark = Tmark(string(Tmark.action) ~= 'E',:);                       % all E events excluded
    Tmark = Tmark(string(Tmark.action) ~= 'R',:);                        % all gaster grooming events excluded
    if height(Tmark)>1 % record the ant only if at least one normal event (not the stroke)
        for i=1:height(Tmark)
            %% work type: A (allo), S (self)
            if string(Tmark.action(i)) == 'G' % allo
                [tf, idx] = ismember(Tmark.R(i), mark); % idx is the index of the receiving ant in mark vector
                loadallo{track_in+idx} = [loadallo{track_in+idx}; [Tmark.frame(i), Tmark.duration(i)]];
                idallofrom{track_in+idmark} = [idallofrom{track_in+idmark}; [track_in+idmark,track_in+idx,Tmark.frame(i),Tmark.duration(i)]];
                idalloto{track_in+idx} = [idalloto{track_in+idx}; [track_in+idmark,track_in+idx,Tmark.frame(i),Tmark.duration(i)]];
            else
                loadself{track_in+idmark} = [loadself{track_in+idmark}; [Tmark.frame(i), Tmark.duration(i)]];
                idself{track_in+idmark} = [idself{track_in+idmark}; [track_in+idmark,track_in+idmark,Tmark.frame(i),Tmark.duration(i)]];
            end
        end
    end
end


for idmark = 1:6
    if length(loadallo{track_in+idmark})>1
        loadallo{track_in+idmark} = sortrows(loadallo{track_in+idmark},1);
    end
    if length(idallofrom{track_in+idmark})>1
        idallofrom{track_in+idmark} = sortrows(idallofrom{track_in+idmark},3);
    end
    if length(idalloto{track_in+idmark})>1
        idalloto{track_in+idmark} = sortrows(idalloto{track_in+idmark},3);
    end
    if length(loadself{track_in+idmark})>1
        loadself{track_in+idmark} = sortrows(loadself{track_in+idmark});
    end
    if length(idself{track_in+idmark})>1
        idself{track_in+idmark} = sortrows(idself{track_in+idmark},3);
    end
end
