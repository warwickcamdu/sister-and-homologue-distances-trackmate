clear all
%% User inputs
% Define path to spots and tracks files
spotsFile='path/to/file';
tracksFile='path/to/file';
% Define sister pairs and homologue groups
sister_pairs = {{0,4}, {5,7}, {10,13}, {14,15}, {20}, {22}};
homologes = {{0,4,5,7}, {10,13,14,15}, {20,22}};

%% Analyse, plot and save data
% Load data
spots = readtable(spotsFile);
tracks = readtable(tracksFile);
%Get base filename
[~, spotsBase, ~]  = fileparts(spotsFile);
spotsBase = erase(spotsBase, '_spots');
% Build output filenames
meanCSV    = sprintf('mean_speeds_%s.csv', spotsBase);
sistersCSV    = sprintf('distance_between_sisters_%s.csv', spotsBase);
homologuesCSV = sprintf('distance_between_homologues_%s.csv', spotsBase);


%% Plot mean speed per track
mean_speeds_table = table(tracks.TRACK_ID, tracks.TRACK_MEAN_SPEED, tracks.TRACK_STD_SPEED, ...
    'VariableNames', {'Track_IDs', 'Mean_Speed', 'Std_Speed'});
writetable(mean_speeds_table,meanCSV);

figure;
errorbar(mean_speeds_table.Track_IDs, mean_speeds_table.Mean_Speed, mean_speeds_table.Std_Speed, 'o');
xlabel('Track ID');
ylabel('Mean Speed');
title('Mean Speed with standard deviation for each Track ID');
grid on;


%% Create subplots for distances
%Find which contain sister track pairs (rather than single tracks)
pairs = {};
for k = 1:length(sister_pairs)
    if length(sister_pairs{k}) == 2
        pairs{end+1} = sister_pairs{k};
    end
end
num_sister_plots = length(pairs);

% Count homologue plots (4-track groups with two sister pairs, and 2-track groups)
num_homologue_plots = length(homologes);

% Set up table
distances_table = table();
% Convert each cell to a string label
sisterNames = cellfun(@(x) sprintf('sisters%dand%d', x{:}), sister_pairs, 'UniformOutput', false);
framesNames = cellfun(@(x) sprintf('frames%dand%d', x{:}), sister_pairs, 'UniformOutput', false);
% Special case: handle single-element cells
for i = 1:numel(sister_pairs)
    if numel(sister_pairs{i}) == 1
        sisterNames{i} = sprintf('sisters%d', sister_pairs{i}{1});
        framesNames{i} = sprintf('frames%d', sister_pairs{i}{1});
    end
end
% Create figures for sisters
figsisters = figure('Name', 'Distances between sisters');
axsisters = gobjects(num_sister_plots,1);
for i = 1:num_sister_plots
    axsisters(i) = subplot(num_sister_plots, 1, i, 'Parent', figsisters);
end
% Plot and save sisters distances
for i = 1:num_sister_plots
    pair = pairs{i};
    [dist, frames, ~] = compute_track_pair_metrics(spots, pair);
    plot_distance(frames, dist, pair, 'Distance between sisters', axsisters(i));
    distances_table.(framesNames{i}) = frames;
    distances_table.(sisterNames{i}) = dist;
end
writetable(distances_table,sistersCSV)

% Create homologue distances table
homologue_table = table();
homologueNames = cellfun(@(x) sprintf('homologues%s', strjoin(string(x), '_')), homologes, 'UniformOutput', false);
framesHomologueNames = cellfun(@(x) sprintf('frames%s', strjoin(string(x), '_')), homologes, 'UniformOutput', false);

% Create figures for homologues
figHomologues = figure('Name', 'Distances between Homologues');
axHomologues = gobjects(num_homologue_plots,1);
for i = 1:num_homologue_plots
    axHomologues(i) = subplot(num_homologue_plots, 1, i, 'Parent', figHomologues);
end

% Plot homologues distances and fill table
hom_idx = 1; % Index for subplot axes
for j = 1:length(homologes)
    if length(homologes{j}) == 4  % Two sister pairs
        dp1 = sister_pairs{2*j-1};
        dp2 = sister_pairs{2*j};
        if length(dp1) == 2 && length(dp2) == 2
            [~, frames1, meanPos1] = compute_track_pair_metrics(spots, dp1);
            [~, frames2, meanPos2] = compute_track_pair_metrics(spots, dp2);

            [commonFrames, distances_between_means] = ...
                compute_mean_to_mean_distance(meanPos1, frames1, meanPos2, frames2);

            plot_distance(commonFrames, distances_between_means, ...
                {dp1{1}, dp1{2}, dp2{1}, dp2{2}}, 'Distance between Homologues', axHomologues(hom_idx));

            % Add to table
            homologue_table.(framesHomologueNames{j}) = commonFrames;
            homologue_table.(homologueNames{j}) = distances_between_means;

            hom_idx = hom_idx + 1;
        end
    elseif length(homologes{j}) == 2  % Simple homolog pair
        [dist, frames, ~] = compute_track_pair_metrics(spots, homologes{j});
        plot_distance(frames, dist, homologes{j}, 'Distance between Homologues', axHomologues(hom_idx));

        % Add to table
        homologue_table.(framesHomologueNames{j}) = frames;
        homologue_table.(homologueNames{j}) = dist;

        hom_idx = hom_idx + 1;
    end
end

% Save table to CSV
writetable(homologue_table, homologuesCSV);

%% Functions
function [distances, commonFrames, meanPositions] = compute_track_pair_metrics(spots, track_ids)
    %Compute distances and mean positions between two tracks
    trackA = spots(spots.TRACK_ID == track_ids{1}, :);
    trackB = spots(spots.TRACK_ID == track_ids{2}, :);

    framesA = trackA.FRAME;
    framesB = trackB.FRAME;
    commonFrames = intersect(framesA, framesB);

    distances = zeros(length(commonFrames), 1);
    meanPositions = zeros(length(commonFrames), 3);

    for i = 1:length(commonFrames)
        frame = commonFrames(i);
        spotA = trackA(trackA.FRAME == frame, :);
        spotB = trackB(trackB.FRAME == frame, :);

        dx = spotA.POSITION_X - spotB.POSITION_X;
        dy = spotA.POSITION_Y - spotB.POSITION_Y;
        dz = spotA.POSITION_Z - spotB.POSITION_Z;
        distances(i) = sqrt(dx^2 + dy^2 + dz^2);

        meanPositions(i, :) = mean([spotA{:, {'POSITION_X','POSITION_Y','POSITION_Z'}};
                                    spotB{:, {'POSITION_X','POSITION_Y','POSITION_Z'}}], 1);
    end
end

function [commonFrames, distances] = compute_mean_to_mean_distance(meanPos1, frames1, meanPos2, frames2)
    % computer the distances between two tracks (or means of two pairs of
    % tracks
    commonFrames = intersect(frames1, frames2);
    distances = zeros(length(commonFrames), 1);

    for i = 1:length(commonFrames)
        f = commonFrames(i);
        idx1 = find(frames1 == f);
        idx2 = find(frames2 == f);

        pos1 = meanPos1(idx1, :);
        pos2 = meanPos2(idx2, :);

        d = pos1 - pos2;
        distances(i) = sqrt(sum(d.^2));
    end
end

function plot_distance(frames, distances, track_ids, label, ax)
    %create distance plot
    plot(ax, frames, distances, '-o');
    xlabel(ax, 'Frame');
    ylabel(ax, 'Distance');
    
    if length(track_ids) == 4
        title(ax, sprintf('%s: Mean of tracks[%d,%d] and Mean of tracks [%d,%d]', ...
            label, track_ids{1}, track_ids{2}, track_ids{3}, track_ids{4}));
    else
        title(ax, sprintf('%s: Track %d and Track %d', ...
            label, track_ids{1}, track_ids{2}));
    end
    grid(ax, 'on');
end
