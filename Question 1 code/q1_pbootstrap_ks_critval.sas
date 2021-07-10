/*GENERATING B SAMPLES FROM H_0 MIXTURE DENSITY TO COMPUTE*/
/*PARAMETRIC BOOTSTRAP CIRITICAL VALUE*/

proc iml;


*K-means subroutine to get inital IM values;

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
			len = ncol(cluster_j); 
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


*_________________________________________________________________________________________________________________________;

n = 200;
pis_h0 = {0.3, 0.25, 0.45};
mus_h0 = {140, 150, 200};
sigmas_h0 = {3, 4, 20};

pis_h1 = {0.2730945, 0.3103271, 0.4165785} ;
mus_h1 = {137.56819, 153.29601, 203.34284};
sigmas_h1 = {3.0408476, 6.5827824, 15.951026};





ncomps = 3; *number of components;
B = 10000; *number of bootstrap KS test statistics we want;

KS_tss = J(B, 1, .);

do bs_iter = 1 to B;

*determining delta for be_iter_i;
delta = {1, 2, 3};
deltas = sample(delta, n, "Replace", pis_h0)`;


x = J(n, 1, .); 
do i = 1 to n;
	delta_i = deltas[i];
	k = delta_i; *while not necessary, this is just to excplicitly show that the value of delta_i governs the component (k) from which we are going to generate an observation.;
	mu_i = mus_h0[k]; *specifiy mean of component k for observation i;
	sigma_i = sigmas_h0[k]; *specifiy standard deviation of component k for observation i;
	x[i] = rand("Normal", mu_i, sigma_i); *the observation generated;
end;

*Now, use this generated sample x, and compute a three component EM density;

points = x;
n = nrow(points);
d = ncol(points); *number of variables in dataset;
k = ncomps; *number of clusters we want;
random = 'obs'; *random inital points selected from observations;
conv = 0; count = 0; pr = 0; *conv must be 0, and will become 1 if the centroids converge. Then, the while loop will stop.;
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


/*INITIAL VALUES FROM K-MEANS*/;

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

mus_em = mus;
sigmas_em = sigmas;
pis_em = pis;


cdf_domain = do(min(x)-20, max(x)+20, 0.1)`;


*cdf of the alternative density in this hypothesis test;
Fx_h0 = 
pis_h0[1]#cdf('normal', cdf_domain, mus_h0[1], sigmas_h0[1]) + 
pis_h0[2]#cdf('normal', cdf_domain, mus_h0[2], sigmas_h0[2]) + 
pis_h0[3]#cdf('normal', cdf_domain, mus_h0[3], sigmas_h0[3]);


Fx_em =
pis_em[1]#cdf('normal', cdf_domain, mus_em[1], sigmas_em[1]) + 
pis_em[2]#cdf('normal', cdf_domain, mus_em[2], sigmas_em[2]) + 
pis_em[3]#cdf('normal', cdf_domain, mus_em[3], sigmas_em[3]);


KS = max(abs(Fx_h0 - Fx_em));

KS_tss[bs_iter] = KS;

end;

*print KS_tss;

/*Calculating the KS test statistic*/

*cdf of the null density in this hypothesis test;
Fx_h0 = 
pis_h0[1]#cdf('normal', cdf_domain, mus_h0[1], sigmas_h0[1]) + 
pis_h0[2]#cdf('normal', cdf_domain, mus_h0[2], sigmas_h0[2]) + 
pis_h0[3]#cdf('normal', cdf_domain, mus_h0[3], sigmas_h0[3]);

*cdf of the alternative density in this hypothesis test;
Fx_h1 = 
pis_h1[1]#cdf('normal', cdf_domain, mus_h1[1], sigmas_h1[1]) + 
pis_h1[2]#cdf('normal', cdf_domain, mus_h1[2], sigmas_h1[2]) + 
pis_h1[3]#cdf('normal', cdf_domain, mus_h1[3], sigmas_h1[3]);


KS_tt = max(abs(Fx_h0 - Fx_h1));
print KS_tt;

/*95% Bootstrap CI for testing a distribution against question 1 H0 density*/
len = nrow(KS_tss);
alpha = 0.1;

call sort(KS_tss);
print KS_tss;

upper_indx = floor((1-(alpha/2))#len);
lower_indx =floor((alpha/2)#len);
ks_rows = nrow(KS_tss);
print ks_rows lower_indx upper_indx;

upper = KS_tss[upper_indx];
lower = KS_tss[lower_indx];

*bootstrap_pvalue;


print "(1-alpha)% CI for KS H0 density critical values" lower upper;











