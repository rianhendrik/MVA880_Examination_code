ods html close;
ods html;

proc iml;
use sasuser.exam_scores;
read all into rd; *rd - reduced digits;
use sasuser.exam_digits;
read all into digits;

true_digits = digits[,1];
df = digits[,2:ncol(digits)]; *df - digits full;

ods listing gpath = "C:\Users\rianh\Dropbox\MVA880\EXAM\LaTeX document\Exam document\figures"
image_dpi=300;


start SetDiag(A, v); *A is the matrix, and v is the value you would like the diagonals to assume;
   diagIdx = do(1,nrow(A)*ncol(A), ncol(A)+1);
   A[diagIdx] = v;             /* set diagonal elements */	
finish;

start kmedoids_clus;

	n = nrow(points);
	d = ncol(points);
	conv = 0;
	count = 0;
	obs = 1:n; *sample the index of the observations for the initial values;
	indices = sample(obs, k, 'WOR');
	inits = points[indices,]; *the inital values, as sampled from the data;

*Assigning each observation to a medoid;
	Fs = 0; R2s = 0;
	iter =1;
	do iter = 1 to 100 while (conv = 0);
		count = count + 1;
		clust_data = inits // points;
		clust_dist = distance(clust_data);
		int = clust_dist[, 1:k]; *int - interest, the part of the distance matrix we are interested in;
		resp_obs = J(nrow(int)-k, 1, .); *subtract k, we are only interested in the datapoints, not the medoids too;
		counter=0;
		do ii = k+1 to nrow(int); 
			counter = counter+1;
			focus = int[ii, ]; resp_obs[counter] = focus[>:<]; *resp - responsibilities; 
		end;

		new_medoids = J(k, d, .); wcd_k = J(k, 1, .);

		do j = 1 to k;
			cluster_j = loc(resp_obs = j);
			points_j = points[cluster_j,]; *values in this cluster;
			len = nrow(points_j); *size of this cluster;
			dist_j = distance(inits[j,]//points_j); *this distance matrix is only used to find the location of the current medoid of the jth cluster;
			call SetDiag(dist_j, 99999);
			current_medoid = loc(dist_j[,1] = 0)-1;
			*Searching for a better mode, with a lower cost using dist_j2;
			dist_j2 = distance(points_j);
			t = inits[j,]; tt = points_j[current_medoid,];
			current_cost = ssq(dist_j2[,current_medoid]);
			other_costs = J(nrow(dist_j2), 1, .);
			do i = 1 to nrow(dist_j2);
				other_costs[i] = ssq(dist_j2[,i]);
			end;
			min_cost_obs = other_costs[>:<]; *the observation which, if it were the medoid, would have the lower cost;
			new_medoids[j,] = points_j[min_cost_obs,];
			wcd_k[j] = other_costs[min_cost_obs];
		end;

	
		
	*Calculating clustering metrics, and checking convergence;
			
		overall_mean = J(1, d, .); *here d denotes dimensions of the dataset;
		do i = 1 to d;
			overall_mean[,i] = points[+,i]/n;
		end;

	*Total deviation;
		td_data = overall_mean // points;
		td_dists = distance(td_data)[, 1]; *We are only interested in the first column, the distance between the centroid and all the individual obervations;
		T = td_dists`*td_dists; *This is the sum of squared distances of each observation from the overall mean; 
		W = sum(wcd_k); B = T - W; tt = W + B;
		c = k;
		r2 = B/T; F = (B/(k - 1))/(W/(n - k)); 
		R2s = R2s//r2; Fs = Fs//F;

		if new_medoids = inits then conv = 1; 
								 else conv = 0;
		if pr = 1 then print iter inits new_medoids ,, W B T tt ,, r2 F;
		inits = new_medoids;
	end;
	
finish;

points = rd;
n = nrow(points);
d = ncol(points); *number of variables in dataset;
k = 10; *number of clusters we want - there are 10 digits;
random = 'obs'; *random inital points selected from observations; *K-means, random points. K-medoids - random from observations;
conv = 0; count = 0; pr = 1; *conv must be 0, and will become 1 if the centroids converge. Then, the while loop will stop.;
call kmedoids_clus;

*Determining the observations in each cluster;

print resp_obs;

i = 2;
t = df[loc(resp_obs = i),];
avrg_gs = t[:,]; *average gray scale per pixel for all image in cluster 1;
clust_i_average = shape(avrg_gs, 16, 16);
title "Heatmap of average pixel grayscale values for digits in cluster 1";
call HeatmapCont(clust_i_average);

