function[transitions,times] = getstats_0D(ant,events,durations);

traces = length(events);            % number of ant traces
transitions = zeros(8,1);           % # transitions X->S, X->A, S->X, S->A, A->X, A->S, S->S, A->A
times = zeros(3,1);                 % times X, S, A

from = {'X','X','S','S','A','A','S','A'};
to = {'S','A','X','A','X','S','S','A'};
states = {'X','S','A'};

%% TRANSITIONS and TIMES 
i = ant;
if length(events{i}) > 1     % only for the relevant ants
    overlaps = find(durations{ant}(2:end,1)-1-durations{ant}(1:end-1,2));    
    if length(overlaps) ~= 0
        %fprintf(1,'\nAnt = %d, Overlaps: %d',ant)
        %fprintf(1,'\nSum times = %1.2f, end of duration = %1.2f',sum(durations{ant}(:,3)),durations{ant}(end,2))
        durations{ant}(overlaps-1:overlaps+1,:);
        events{ant}(overlaps-1:overlaps+1,:);
        durations{ant}(overlaps+1,:) = [];
        events{ant}(overlaps+1) = [];
    end
    LE = length(events{i});         % number of states in the current ant trace
    TE = durations{i}(:,1);         % starting time of the event
    
    
    if length(durations{i}) > 0
        T0 = ceil(durations{i}(1,1));     % starting frame (0)
        T1 = ceil(durations{i}(end,2));   % ending frame (~30 000)
        Tvec = [T0:T1];             % frame vector
        Svec = [];                  % vector of states
        Svec = repelem(events{i},ceil(durations{i}(:,3)))';
        
    end
    
    %% information on each transition
    Sfrom = events{i}(1:LE-1);
    Sto = events{i}(2:LE);
    for react = 1:8
        nr = length(find(Sfrom == from(react) & Sto == to(react))); %% numbers of different transitions
        transitions(react) = transitions(react) + nr;
    end
    
    
    
    %% count times at each state    
    Sfrom = Svec(1:end);
    for is = 1:3
        ns = length(find(Sfrom == states(is))); %% indices of different transitions
        times(is) = times(is) + ns;
    end
end


