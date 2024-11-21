%Vikram Vijayakumar (02068559)
%MTH 565 Project 6_1
%https://www.mathworks.com/help/matlab/ref/loglog.html
%I had reviewed the above github link to obtain examples for log functions

t = 500;  %Number of time steps, adjust for 100 and 500
initial_nodes = 3;  %Initial number of nodes

%Generate the graph G1(t)
G = zeros(initial_nodes); %Initial graph (3 nodes, 2 edges)
G(1, 2) = 1; G(2, 1) = 1;  %Edge 1-2
G(2, 3) = 1; G(3, 2) = 1;  %Edge 2-3
degrees = sum(G, 2); %Initial degrees

for step = 1:t
    % Add a new node
    new_node = size(G, 1) + 1; % Add a new node
    G(new_node, :) = 0;  % Initialize new row
    G(:, new_node) = 0;  % Initialize new column
    
    %Connect to an existing vertex based on preferential attachment
    probabilities = degrees / sum(degrees); 
    target_node = randsample(1:length(degrees), 1, true, probabilities);
    G(new_node, target_node) = 1;
    G(target_node, new_node) = 1;
        
%add an extra edge between two randomly chosen distinct vertices
if rand < 1 %Adjust probability as needed
    eligible_nodes = find(degrees > 0); % Only nodes with non-zero degree
    if length(eligible_nodes) > 1 % Ensure there are at least two eligible nodes
        %Select two distinct vertices
        [v1, v2] = deal(0, 0);
        while v1 == v2 || G(v1, v2) == 1
            v1 = randsample(eligible_nodes, 1, true, degrees(eligible_nodes));
            v2 = randsample(eligible_nodes, 1, true, degrees(eligible_nodes));
        end
        %Add edge between the selected vertices
        G(v1, v2) = 1;
        G(v2, v1) = 1;
    end
end
        
%Update degrees
degrees = sum(G, 2);
end

% Plot the graph
figure;
graph_obj = graph(G); % Convert adjacency matrix to a graph object
plot(graph_obj, 'Layout', 'force'); % Use a force-directed layout for visualization
title('Graph Visualization of G1(t)');


%Calculate diameter
diameter = calculate_diameter(G);
fprintf('Diameter of G1(%d): %d\n', t, diameter);

%Calculate clustering coefficients
clustering_coefficients = calculate_clustering_coefficient(G);
avg_clustering = mean(clustering_coefficients(clustering_coefficients > 0));
fprintf('Average Clustering Coefficient of G1(%d): %.4f\n', t, avg_clustering);

%Plot histogram of vertex degrees
degrees = sum(G, 2); % Degree of each node
figure;
histogram(degrees, 'Normalization', 'probability');
title('Histogram of vertex degrees');
xlabel('Degree');
ylabel('Probability');

%Check for power-law distribuition
deg_counts = tabulate(degrees);
large_degrees = deg_counts(deg_counts(:,1) > 2, :); %consider only larger degrees
log_x = log(large_degrees(:,1));
log_y = log(large_degrees(:,2));
coeffs = polyfit(log_x, log_y, 1);

% Plot the power-law fit
figure;
plot(log_x, log_y, 'o');
hold on;
plot(log_x, polyval(coeffs, log_x), '-r');
title('Power-Law for vertex degree');
xlabel('Degree');
ylabel('(Frequency');

%Function to calculate diameter
function d = calculate_diameter(G)
    graph_obj = graph(G);
    D = distances(graph_obj); % Shortest path distances
    D(D == Inf) = 0; % Exclude disconnected pairs
    d = max(D(:)); % Diameter is the longest shortest path
end

%Function to calculate clustering coefficient
function C = calculate_clustering_coefficient(G)
    n = size(G, 1);
    local_cluster = zeros(1, n);

    for i = 1:n
        neighbors = find(G(i, :));
        k = length(neighbors);
        if k > 1
            subG = G(neighbors, neighbors);
            num_links = sum(subG(:)) / 2;
            local_cluster(i) = (2 * num_links) / (k * (k - 1));
        else
            local_cluster(i) = 0;
        end
    end
    C = mean(local_cluster);
end