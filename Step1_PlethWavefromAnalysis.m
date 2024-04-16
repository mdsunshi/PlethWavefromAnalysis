clc
clear
close all hidden

%% files
files = ; %generate list of files

%%
for f = 1:numel(files)
    cd(file_path)
    fname = files(f).name;
    load(fname);% will need to edit to fit the file type
    data = ; % assign trace to variable data

    %%
    %% Assuming the data is in a single column, if not, adjust accordingly
    time = (0:length(data)-1) / 1000; % Time vector in seconds

    %% Design a bandpass filter for 0.1-10 Hz
    fs_original = 1000; % Original sampling rate
    fs_target = 90; % Target sampling rate
    nyquist = fs_original / 2;
    filter_high = 2 / nyquist;
    filter_low = 10 / nyquist;
    [b, a] = butter(2, [filter_high, filter_low]);

    %% Apply the bandpass filter
    filtered_data = filtfilt(b, a, data);

    %% Resample data to 90 samples per second
    resampled_data = resample(filtered_data, fs_target, fs_original);

    %% find breaths
    MPP = 0.5;
    MPH = 0.5;
    MPD = 0.1*fs_target;
    [pks,locs] = findpeaks(resampled_data,...
        'MinPeakProminence',MPP,...
        'MinPeakHeight',MPH,...
        'MinPeakDistance',MPD);

    valley = zeros(size(locs));
    beforepeak = 0.5*fs_target;
    for i = 1:numel(locs)
        start = floor(locs(i)-beforepeak);
        stop = (locs(i));
        range = start:stop;

        if (range(1) <= 0) || (range(end) > numel(resampled_data))
            locs(i) = 0;
            continue
        else
        end

        TF = islocalmin(resampled_data(range));
        pos = find(TF);
        if ~isempty(pos)
            valley(i) = range(pos(end));
        else
            locs(i) = 0;
            continue
        end
    end
    valley = valley(find(valley));
    locs = locs(find(locs));

    % figure
    % hold on
    % plot(resampled_data)
    % scatter(locs,resampled_data(locs),10,'rv','filled')
    % scatter(valley,resampled_data(valley),10,'g^','filled')

    BreathMat = zeros(numel(valley)-1,880);
    for i = 1:numel(valley)-1
        breath_bit = resampled_data(valley(i):valley(i+1));
        breath_bit_cycleNorm = resample(breath_bit, 900, numel(breath_bit));
        breath_bit_cycleNorm = breath_bit_cycleNorm(1:880);
        
        % plot(valley(i):valley(i+1),breath_bit+0.2)
        % plot(linspace(valley(i),valley(i+1),numel(breath_bit_cycleNorm)),breath_bit_cycleNorm+0.2)
        

        % BreathCell{i} = breath_bit;
        BreathMat(i,:) = breath_bit_cycleNorm;
    end

    %% save data matrix for each file
    cd(save_path)
    save(['BreathMat_' fname],'BreathMat','valley','fs_target')

end


