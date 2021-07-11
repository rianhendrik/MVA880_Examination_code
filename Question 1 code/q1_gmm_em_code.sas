ods html close;
ods html;

proc iml;
use sasuser.exam_dist;
read all into x;

ncomps = 3; *the number of components we believe there to be;

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
k = ncomps; *number of clusters we want;
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

/*INITIAL VALUES*/;
*as stated in my document, these intial values were obtained from my k-means algorithm, also available in my submission document;

mus = means;
sigmas = stds;
pis = pis;

n = nrow(x);
tol = 0.01; *tollerance for convergence;
z = 2; *to start of the while loop;
liks = {0, 1}; *to start of the while loop's condition;
iter = 0; *to keep track of iterations before convergence.;

do while (abs(liks[z] - liks[z - 1]) >= tol);
	
	iter = iter + 1;
	
	*calculating the loglik to know when iterations should stop;
	


	/*E Step*/
	comps = J(n, ncomps, .);
	do k = 1 to ncomps;
		comps[,k] = pis[k,]*pdf('normal', x, mus[k,], sigmas[k,]);
	end;

	comp_sum = comps[,+]; *for each row, sum accross each column - the ith row of this, gives us the denominator of gamma_ik for the ith obervation;
	gammas = comps/comp_sum; *calculating our gamma_ik's
	
	/*M Step*/;

	*effective sample size;
	nks = gammas[+,]; *our effective sample sizes;

	*pi's;
	pis_n = (nks/n)`; *our new pi's, after iteration z - the k'th effective sample size divided by the total number of observation = k'th pi.;

	*mu's;
	mus_n = (((gammas#x)[+,])/nks)`;


	*simga's;
	*centering each x, and squaring each centred x k times using the kth mean, to compute k standard deviations;
	xs = J(n, ncomps, .);
	do iter1 = 1 to ncomps;
		xs[,iter1] = (x - mus_n[iter1])##2;
	end;
	sigmas_n = sqrt(t(vecdiag(gammas`*xs))/nks)`; 

	z = z + 1;

	*if pr = 1 then print iter ,, pis_n mus_n sigmas_n;
	*print iter ,, pis_n mus_n sigmas_n;


	
	*updating our parameter estimates after iteration z;
	pis = pis_n; 
	mus = mus_n;
	sigmas = sigmas_n;

	lik = J(n, ncomps, .);
	do iter0 = 1 to ncomps;
		lik[,iter0] = pis[iter0,]#pdf('normal', x, mus[iter0,], sigmas[iter0,]);
	end;

	loglik = sum(log(lik[,+]));
			
	liks = liks//loglik;

end;

print liks;

BIC = ncomps#3*log(n) - 2#liks[nrow(liks)];
AIC = ncomps#3#2 - 2#liks[nrow(liks)];
print BIC,,AIC;


*Determining the hard clustering solution;

clusts = J(n, 1, .);
do i = 1 to n;
	clusts[i] = loc(gammas[i,] = max(gammas[i,]));
end;

*Sorting cluster numbers, so that cluster with samllest mean is #1, etc;
map = {1 2,	
	   2 3,
       3 1};
	
sorted_clusts = J(n, 1, .);
do i = 1 to n;
	clusti = clusts[i];
	do j = 1 to 3;
		if clusti = map[j,1] then sorted_clusts[i] = map[j,2];
	end;
end;

print clusts sorted_clusts x;

hard_clus_solution = x||sorted_clusts;

create gmm_hard_clus_solution from hard_clus_solution;
append from hard_clus_solution;

quit;
