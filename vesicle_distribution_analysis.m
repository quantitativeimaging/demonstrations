% Plot vesicle nearest-neighbour distance histogram
%  EJR 2017 cc-by
%
% Also compare observed and simulated (from uniform random dist) number
% of nearest neighbours within a critical radius of interest.
% If an adhesion exists, the observed number may be > simulated
%
% Notes: 
%   Dimensions = pixel widths
%   Analysis depends on particle detection working correctly to produce
%  x- and y- positions. Outside scope of this script. 
%
% Sample data:
% D:\EJR_OneDrive\OneDrive - University Of Cambridge\Projects\2017_vesicle_distribution\EGTA no aSyn

[filename_in,pathname_in] = uigetfile('*.xls','Select the XLS particle location file');

dist_threshold = 10; % units: pixel widths

data_xls = xlsread([pathname_in, filename_in]);

x = data_xls(21:end, 3);
y = data_xls(21:end, 4);

interparticle_distance_matrix = squareform( pdist([x,y],'euclidean') );
interparticle_distance_matrix((eye(length(x)))==1) = inf; % ignore self

nearest_neighbour_distances = min(interparticle_distance_matrix, [], 2);

% At this point we could just plot the nearest neighbour histograms for 
% each treatment and compare them
figure(1)
hist(nearest_neighbour_distances, [0:5:100])
xlabel('nearest neighbour distance, pixel widths', 'fontSize', 14)
ylabel('number', 'fontSize', 14)
set(gca, 'fontSize', 14)
set(gcf, 'color', 'white')
xlim([0 100])

Number_sub_crit_distance = sum(nearest_neighbour_distances < dist_threshold);

% Simulate expected number within this distance (try 20 simulations):
list_sim_number_sub_crit_distance = zeros(20,1);
for lp = 1:20
    simX = rand(size(x))*1000;
    simY = rand(size(y))*1000;
    
    sim_interparticle_distance_matrix = squareform( pdist([simX,simY],'euclidean') );
    sim_interparticle_distance_matrix((eye(length(x)))==1) = inf;

    sim_nearest_neighbour_distances = min(sim_interparticle_distance_matrix, [], 2);
    sim_number_sub_crit_distance = sum(sim_nearest_neighbour_distances < dist_threshold);
    list_sim_number_sub_crit_distance(lp) = sim_number_sub_crit_distance;
end

Expt = Number_sub_crit_distance
Sim_mean  = mean(list_sim_number_sub_crit_distance)
Sim_std   = std(list_sim_number_sub_crit_distance)

