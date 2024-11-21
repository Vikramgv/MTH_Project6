%Vikram Vijayakumar (02068559)
%MTH 565 Project 6_2
%https://www.mathworks.com/help/matlab/ref/loglog.html
%I had reviewed the above github link to obtain examples for log functions

t = 100; %Number of time steps, adjust for 100 and 500
p = 0.5; %Probability of adding a new vertex
initial_nodes = 3; % Initial number of nodes
max_nodes = initial_nodes + t; % Maximum possible nodes

% Preallocate a sparse matrix
G = sparse(max_nodes, max_nodes);

% Initialize the initial graph (3 nodes, 2 edges)
G(1, 2) = 1; G(2, 1) = 1; % Edge 1-2
G(2, 3) = 1; G(3, 2) = 1; % Edge 2-3

% Initial degrees and total degree
degrees = sum(G, 2);
total_degree = sum(degrees);

% Generate Gp(t)
for step = 1:t
    if rand <= p
% Add a new vertex
new_node = initial_nodes + step;

% Ensure probabilities are correctly computed
if total_degree == 0
    probabilities = ones(initial_nodes + step - 1, 1) / (initial_nodes + step - 1);
else
    probabilities = degrees(1:initial_nodes + step - 1) / total_degree;
end

% Preferential attachment based on probabilities
target_node = randsample(1:initial_nodes + step - 1, 1, true, probabilities);

% Add the new edge
G(new_node, target_node) = 1;
G(target_node, new_node) = 1;

% Update degrees and total degree
degrees(new_node) = 1;
degrees(target_node) = degrees(target_node) + 1;
total_degree = total_degree + 2;

    end
end

% Trim unused rows and columns
G = G(1:initial_nodes + t, 1:initial_nodes + t);

% Plot graph
graph_obj = graph(G);
plot(graph_obj, 'Layout', 'force');

% Calculate diameter
diameter = calculate_diameter(Gp);
fprintf('Diameter of Gp(%d): %d\n', t, diameter);

% Calculate clustering coefficients
clustering_coefficients = calculate_clustering_coefficient(Gp);
avg_clustering = mean(clustering_coefficients(clustering_coefficients > 0));
fprintf('Average Clustering Coefficient of Gp(%d): %.4f\n', t, avg_clustering);

% Plot histogram of vertex degrees
degrees = sum(Gp, 2); % Degree of each node
figure;
histogram(degrees, 'Normalization', 'probability');
title('Degree Distribution of Gp(t)');
xlabel('Degree');
ylabel('Probability');

% Check for power-law behavior in the tail
deg_counts = tabulate(degrees);
large_degrees = deg_counts(deg_counts(:, 1) > 2, :); % Focus on large degrees
log_x = log(large_degrees(:, 1));
log_y = log(large_degrees(:, 2));

% Fit a linear model to log-log data
coeffs = polyfit(log_x, log_y, 1);
slope = coeffs(1);
intercept = coeffs(2);
fprintf('Power-law fit for Gp(t): y = %.4f * x^(%.4f)\n', exp(intercept), slope);

% Plot the power-law fit
figure;
plot(log_x, log_y, 'o');
hold on;
plot(log_x, polyval(coeffs, log_x), '-r');
title('Power-Law Fit of Degree Distribution for Gp(t)');
xlabel('Degree');
ylabel('Frequency');

% Function to calculate diameter
function d = calculate_diameter(G)
    graph_obj = graph(G);
    D = distances(graph_obj); % Shortest path distances
    D(D == Inf) = 0; % Exclude disconnected pairs
    d = max(D(:)); % Diameter is the longest shortest path
end

% Function to calculate clustering coefficient
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