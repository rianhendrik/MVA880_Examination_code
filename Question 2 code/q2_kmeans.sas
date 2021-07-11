ods html close;
ods html;

proc iml;
use sasuser.exam_scores;
read all into rd; *rd - reduced digits;
use sasuser.exam_digits;
read all into digits;

ods listing gpath = "C:\Users\rianh\Dropbox\MVA880\EXAM\LaTeX document\Exam document\figures"
image_dpi=300;

true_digits = digits[,1];
d = digits[,2:ncol(digits)];
print true_digits;

/*r = nrow(rd);*/
/*c = ncol(rd);*/
/*print r c;*/

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
		counter = 0;
		do ii = k+1 to nrow(int); *start at k, we are not clustering the k centroids.;
			counter = counter+1;
			focus = int[ii, ]; resp[counter] = focus[>:<]; *resp - responsibilities;
		end;
		
		print resp;
		

		new_centroids = J(k, d, .); wcd_k = J(k, 1, .);

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

*print test;

points = rd;
n = nrow(points);
d = ncol(points); *number of variables in dataset;
k = 10; *number of clusters we want - there are 10 digits;
random = 'obs'; *random inital points selected from observations; *K-means, random points. K-medoids - random from observations;
*inits = {95 95, 100 100, 105 105}; *predetermined inital values 2D;
conv = 0; count = 0; pr = 1; *conv must be 0, and will become 1 if the centroids converge. Then, the while loop will stop.;
call kmeans_clus;
print f r2;

*Determining the observations in each cluster;

print resp_obs;

i = 10;
t = d[loc(resp_obs = i),];
avrg_gs = t[:,]; *average gray scale per pixel for all image in cluster 1;
clust_i_average = shape(avrg_gs, 16, 16);
title "Heatmap of average pixel grayscale values for digits in cluster 1";
call HeatmapCont(clust_i_average);

*note that these allocations will change the next time the k-means algorithm is run;
kmeans_class = resp_obs;
do i = 1 to n;
	if resp_obs[i] = 6 then kmeans_class[i] = 0;
	if resp_obs[i] = 1 then kmeans_class[i] = 1;
	if resp_obs[i] = 10 then kmeans_class[i] = 2;
	if resp_obs[i] = 9 then kmeans_class[i] = 3;
	if resp_obs[i] = 2 then kmeans_class[i] = 4;
	if resp_obs[i] = 4 then kmeans_class[i] = 5;
	if resp_obs[i] = 8 then kmeans_class[i] = 6;
	if resp_obs[i] = 3 then kmeans_class[i] = 7;
	if resp_obs[i] = 7 then kmeans_class[i] = 8;
	if resp_obs[i] = 5 then kmeans_class[i] = 9;
end;

kmeans_res = true_digits[41:n]||kmeans_class[41:n];

res = kmeans_res;

create kmeans from res;
append from res;

quit;


*determining the proportion of data misclassified;
proc iml;
use sasuser.exam_kmeans_class;
read all into res;

truth = res[,1];
kmc = res[,2];

missclass = J(nrow(res), 1, 1);
do i = 1 to nrow(res);
	if truth[i] = kmc[i] then missclass[i] = 0;
end;

prop_miss = sum(missclass)/nrow(missclass);
print prop_miss missclass;
	





