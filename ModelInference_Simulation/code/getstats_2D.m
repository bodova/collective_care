function[transitions,times,histY12] = getstats_2D(ant,var1,var2,events,durations,M,N1,y1,N2,y2,DTA,DTL,ifL);
global idty Lreset

% var1=charACT_P;
% var2=loadSEEN_a;
% M=[0,0];
% N1=NP; y1=yP; N2=NL; y2=yL;
% ifL=[0,1];

traces = length(events);            % number of ant traces
transitions = zeros(8,N1,N2);       % # transitions X->S, X->A, S->X, S->A, A->X, A->S, S->S, A->A
times = zeros(3,N1,N2);             % times X, S, A
histY12 = [];                        % activity percentage

%% y1, y2 - bin centers, get bin edges
d = diff(y1)/2; y1edges = [y1(1)-d(1), y1(1:end-1)+d, 10^6];
d = diff(y2)/2; y2edges = [y2(1)-d(1), y2(1:end-1)+d, 10^6];

from = {'X','X','S','S','A','A','S','A'};
to = {'S','A','X','A','X','S','S','A'};
states = {'X','S','A'};

%fprintf(1,'\nComputing transition statistics...(%d)',traces)
%% TRANSITIONS and TIMES depending on explanatory variable Y1 - activity performed, Y2 - activity received
%for i = 1:traces % number of ant traces
i = ant;
if length(events{i}) > 1     % only for the relevant ants
    
    overlaps = find(durations{ant}(2:end,1)-1-durations{ant}(1:end-1,2));
    if length(overlaps) ~= 0
        if (durations{ant}(overlaps,2)>=durations{ant}(overlaps+1,1)) & (durations{ant}(overlaps,2)<durations{ant}(overlaps+1,2))
            durations{ant}(overlaps,2) = durations{ant}(overlaps+1,1) - 1;
            durations{ant}(overlaps,3) = durations{ant}(overlaps,2) - durations{ant}(overlaps,1) + 1;
        else
            durations{ant}(overlaps+1,:) = [];
            events{ant}(overlaps+1) = [];
        end
    end
    
    % Activity data: charACT or charACTself - which activity matters
    A = var1{i};
    B = var2{i};

    
    [Tvec,Svec,Bvec,Btran] = detailed_states_trace(i,events,durations,B,ifL(2),DTL, Lreset); 
    %% Tvec - vector of frames, Svec - vector of states in each frame
    LE = length(events{i});         % number of states in the current ant trace
    TE = durations{i}(:,1);         % starting time of the event
    
    
    
    if ifL(1) == 0 
        Amean = [0;movmean(A,[DTA-1,0])];
    end
    if ifL(2) == 0 
        Bmean = [0;movmean(B,[DTA-1,0])];
    end
    
    
    %%%%%%% assume only the second variable can be L
    if idty == 2 
        Y1 = Amean(TE(2:LE)+1);
    else
        if TE(LE)+1 < length(Amean)
            Y1 = Amean(TE(2:LE)+1) - M(1);
        else
            Y1 = Amean(TE(2:LE-1)+1) - M(1);
        end
    end
        
        
    if idty == 2 && ifL(2) == 0
        Y2 = Bmean(TE(2:LE)+1);
    elseif ifL(2) == 0
        if TE(LE)+1 < length(Bmean)
            Y2 = Bmean(TE(2:LE)+1) - M(2);
        else
            Y2 = Bmean(TE(2:LE-1)+1) - M(2);
        end
    elseif ifL(2) == 1
        Y2 = Bvec(durations{i}(2:end,1))';
    end
    
    
    
    %% information on each transition: which bin, state from and to
    Sfrom = events{i}(1:LE-1);
    Sto = events{i}(2:LE);
    LY1 = length(Y1); LY2 = length(Y2); LY = min(LY1,LY2); 
    [N,edgeX,edgeY,binX,binY]=histcounts2(Y1(1:LY),Y2(1:LY),y1edges,y2edges); 
    
    %% collect all transitions
    for react = 1:8
        ir = find(Sfrom == from(react) & Sto == to(react)); %% indices of different transitions
        ir = ir(ir<=LY);
        [counts,h] = hist3([binX(ir),binY(ir)],{1:length(y1) 1:length(y2)});
        for ibin = 1:length(y1)
            for jbin = 1:length(y2)
                transitions(react,ibin,jbin) = transitions(react,ibin,jbin) + counts(ibin,jbin);
            end
        end
    end
    
    %% collect all times
    if idty == 2
        Y1 = Amean;
    else
        Y1 = Amean;
    end
    
    if idty == 2 && ifL(2) == 0
        Y2 = Bmean;
    elseif ifL(2) == 0
        Y2 = Bmean - M(2);
    else
        Y2 = Bvec(Tvec(1:end)+1)';
    end

        
    LY1 = length(Y1); LY2 = length(Y2); LY = min([LY1,LY2,length(Svec)]);   
    Sfrom = Svec(1:LY);
    [N,edgeX,edgeY,binX,binY]=histcounts2(Y1(1:LY),Y2(1:LY),y1edges,y2edges);
    
    for is = 1:3
        ir = find(Sfrom == states(is)); %% indices of different states
        [counts,h] = hist3([binX(ir),binY(ir)],{1:length(y1) 1:length(y2)});
        for ibin = 1:length(y1)
            for jbin = 1:length(y2)
                times(is,ibin,jbin) = times(is,ibin,jbin) + counts(ibin,jbin);
            end
        end
    end
    histY12 = [histY12;[Y1(1:LY),Y2(1:LY)]];          % collects (Y1,Y2) joint values
    
    
end



