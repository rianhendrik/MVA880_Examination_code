ods html close;
ods html;

proc iml;
use sasuser.exam_dist;
read all into x;


start kmeans_clus;

	n = nrow(points);
	d = ncol(points);
	conv = 0;
	count = 0;
	obs = 1:n; *sample the index of the observations for the initial values;
	if random = 'obs' then indices = sample(obs, k, 'WOR'); if random = 'obs' then inits = points[indices,];

	overall_mean = J(1, d, .); *here d denotes dimensions of the dataset;
	do i = 1 to d;
		overall_mean[,i] = points[+,i]/n;
	end;

*Total deviation;
	td_data = overall_mean // points;
	td_dists = distance(td_data)[, 1]; *We are only interested in the first column, the distance between the centroid and all the individual obervations;
	T = td_dists`*td_dists; *This is the sum of squared distances of each observation from the overall mean; 

*Assigning each observation to a centroid;
	Fs = 0; R2s = 0;
	do iter = 1 to 100 while (conv = 0);
		count = count + 1;
		clust_data = inits // points;
		clust_dist = distance(clust_data);
		int = clust_dist[, 1:k]; *int - interest, the part of the distance matrix we are interested in;
		resp = J(nrow(int), 1, .);
		do ii = i to nrow(int); *start at k, we are not clustering the k centroids.;
			focus = int[ii, ]; resp[ii] = focus[>:<]; *resp - responsibilities;
		end;
		resp_obs = resp[k + 1:nrow(resp),];

		new_centroids = J(k, d, .); wcd_k = J(k, 1, .);
		*print iter;
		do j = 1 to k;
			cluster_j = loc(resp_obs = j);
			len = ncol(cluster_j); print len;
			points_j = points[cluster_j,];
			wcd_data = inits[j, ] // points_j; wcd_dist = distance(wcd_data)[, 1]; *We are onlt interested in the distane of each cluster centroid to the observations in its cluster;
			wcd_k[j] = wcd_dist`*wcd_dist;
			new_centroids[j, ] = points_j[+,]/nrow(points_j);
		end;
		W = sum(wcd_k); B = T - W; tt = W + B;
		c = k;
		r2 = B/T; F = (B/(k - 1))/(W/(n - k)); 
		R2s = R2s//r2; Fs = Fs//F;

		if new_centroids = inits then conv = 1; 
								 else conv = 0;
		if pr = 1 then print iter inits new_centroids ,, W B T tt ,, r2 F;
		inits = new_centroids;
	end;
finish;

points = x;
n = nrow(points);
d = ncol(points); *number of variables in dataset;
k = 2; *number of clusters we want;
random = 'obs'; *random inital points selected from observations;
conv = 0; count = 0; pr = 1; *conv must be 0, and will become 1 if the centroids converge. Then, the while loop will stop.;
call kmeans_clus;

*Calculating mean and standard deviation and relative size of the clusters - for GMM initial values;
solution = points || resp_obs;

means = J(k, 1, .);
stds = J(k, 1, .);
pis = J(k, 1, .);
do i = 1 to k;
	clust_i = solution[loc(solution[,2] = i),1];
	means[i] = mean(clust_i);
	stds[i] = std(clust_i);
	pis[i] = nrow(clust_i)/n;
end;

print means stds pis;
	
	



