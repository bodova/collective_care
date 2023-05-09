function[Tvec,Svec,Lvec,Ltran] = detailed_states_trace(i,events,durations,A,ifL,DTL,Lreset);
%% Tvec - vector of frames
%% Svec - vector of states in each frame

if length(durations{i}) > 0
    T0 = ceil(durations{i}(1,1));     % starting frame (0)
    T1 = ceil(durations{i}(end,2));   % last ending frame (~30 000)
    T2 = ceil(durations{i}(end-1,2));   % previous ending frame (~30 000)
    Tvec = [T0:T1];             % frame vector
    Svec = [];                  % vector of states
    Lvec = [];
    Ltran = [];
    L = length(events{i});
    
    Svec = repelem(events{i},ceil(durations{i}(:,3)))';
    
    if ifL == 1 % computes Ltran - values of L at the transition times
                % computes Lvec - statistics, which is equal to the last
                % nonzero seen load (not your own!!!) for every time in

        Ltran = 0;
        for j = 2:L-1 % number of recorded events in one trace
            if Lreset == 1 % only considers last DT of data
                jj = durations{i}(j,1);
                if jj-DTL < 1
                    Llast = A(find(A(1:jj),1,'last'));
                elseif max(A(jj-DTL:jj))==0
                    Llast = 0;
                else
                    Llast = A(find(A(1:jj),1,'last'));
                end
            else % considers the whole trace 
                Llast = A(find(A(1:durations{i}(j,1)),1,'last')); % last load seen from the beginning of the experiment
            end
            if length(Llast) == 0
                Llast = 0;
            end
            Ltran = [Ltran; Llast];
            %fprintf(1,'\nj=%d, Time: %d--%d, Lend=%1.2f, (%c,%c), Llast=%1.2f\n',j,max(0,durations{i}(j,1)-DTL),durations{i}(j,1),A(durations{i}(j,1)),events{i}(j-1),events{i}(j),Llast)

        end
        
        
        
        %% we consider load to change during the states 
        Lvec = [];
        for j = 2:T1 % all times
            if Lreset == 1 % only considers last DT of data
                if j-DTL < 1
                    Llast = A(find(A(1:j),1,'last'));
                elseif j >= length(A)
                    Llast = A(find(A,1,'last'));
                elseif max(A(j-DTL:j))==0
                    Llast = 0;
                else
                    Llast = A(find(A(1:j),1,'last'));
                end
            else % considers the whole trace 
                Llast = A(find(A(1:j),1,'last')); % last load seen from the beginning of the experiment
            end
            
            if length(Llast) == 0
                Llast = 0;
            end
            Lvec(j) = Llast;
        end
        
        id1 = (A(1:min(T1+1,length(A)))>0); % all times when allo action is done and the load is positive
        Lvec(id1) = A(id1);
        Lvec = [0,Lvec]; % for the dimensions to be consistent with the effective approach
    end

    
else
    Tvec = [];
    Svec = [];
    Lvec = [];
end

