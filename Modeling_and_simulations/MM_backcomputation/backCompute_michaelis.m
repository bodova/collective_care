% Reads from a file ('POST_windowed_450_groomTot.csv', see Preparation in this 
% repository to reproduce this table) which contains one row per treated
% ant. The first 9 columns describe the ant and measured spores, the next numWindows
% are float with the number of frames of total grooming received (i.e. allogrooming 
% and body selfgrooming) by treated ants up to the start of each window.
% The output has the same format, the first 9 columns being identical as in
% the input, and the next numWindows columns are a time series representing
% the estimated number of body spores at the start of the window, given by the
% Michaelis-Menten equation.

% Example execution 
%[binsload extrapinitial] = backCompute_michaelis(5.22, 50119, 'POST_windowed_450_groomTot_spore_MM.csv')

%%It reads from a file( /MM_backcomputation/'POST_windowed_450_groomTot.csv', provided but see Preparation in this repository to reproduce this table from raw data) which contains one row per treated ant, and includes a windowed time series of total grooming received (i.e. allogrooming and body selfgrooming) by treated ants up to the start of each window.
The output has the same format and the time series contains the estimated number of body spores at the start of the window, given by the Michaelis-Menten equation.
Example execution is given in the header.


function [binsload extrapinitial] = backCompute_michaelis(v, K, outputfile)
    %numWindows=180;

    %% parse the input file
    fid = fopen('./POST_windowed_450_groomTot.csv');

    str_init='%q%f%q%q%q%q%q%q%q';
    str = str_init;
    for i=1:numWindows,
        str = [str '%f'];
    end

    CC = textscan(fid,str,'delimiter',',');
    fclose(fid);

    tau_data = zeros(numel(CC{1}),numWindows);

    for i=1:numel(CC{1}),
        tx = CC{4}(i);
        for q=1:numWindows,
            tau_data(i,q) = CC{9+q}(i);
        end
        if (strcmp(tx{1},'R'))
            fx = CC{8}(i);
            fload(i) = str2double(fx);
        elseif (strcmp(tx{1},'G'))
            fx = CC{9}(i);
            fload(i) = str2double(fx);
        else
            fload(i) = nan;
            tau_data(i,:) = nan;
        end
    end

    %% Compute number of spores
    % back integrate dS/dt = -v S / (S+K), from S_final at T_f
    binsload = zeros(size(tau_data,1),size(tau_data,2)+1);
    binsload(:,end) = fload;
    for q=(numWindows):-1:1,
        binsload(:,q) = binsload(:,q+1) + v * binsload(:,q+1) ./ (binsload(:,q+1)+K) .* tau_data(:,q);
    end

    %% Write output file
    str_init = strrep(str_init,'q','s');
    str_init = strrep(str_init,'%',',%');
    str_init = strrep(str_init,'f','d');
    str_init = str_init(2:end);
    fidout = fopen(outputfile, 'w');
    for i=1:numel(CC{1})
        fprintf(fidout,str_init,CC{1}{i},CC{2}(i),CC{3}{i},CC{4}{i},CC{5}{i},CC{6}{i},CC{7}{i},CC{8}{i},CC{9}{i});
        fprintf(fidout,',%f',binsload(i,:));
        fprintf(fidout,'\n');
    end
    fclose(fidout);
    

    extrapinitial = [];
    % now explicitly solve differential equation solution for the initial
    % load, and compare with backintegration; extrapinitial should be same
    % as binsload(:,1)...
%     tau = sum(tau_data');
%     for q=1:numel(CC),
%         if (isnan(fload(q)))
%             extrapinitial(q) = nan;
%         else
%             extrapinitial(q) = fsolve(@(x)(fload(q) + v * tau(q) + K * log(fload(q) / x) - x), 100000);
%         end
%     end
