/*COMPARING K-NN CLUSTERING SOLUTION TO GMM HARD CLUSTERING SOLUTION*/


/*SPECIFYING PATH FOR PLOTS, AND PLOT QUALITY*/
ods listing gpath = "C:\Users\rianh\Dropbox\MVA880\EXAM\LaTeX document\Exam document\figures"
image_dpi=300;


ods html close;
ods html;

proc iml;
use sasuser.exam_knn_clus_solution; read all into knn; knn =knn[,1:2];
use sasuser.exam_gmm_hard_clus_solution; read all into gmm;
use sasuser.exam_knn_density_est; read all into density;

modes = unique(knn[,2]);
print modes;
clusts = J(nrow(knn), 1, .);
do i = 1 to nrow(knn);
	do j = 1 to 3;
		if knn[i, 2] = modes[j] then clusts[i] = j;
	end;
end;

knn_solution = knn[,1]||clusts; print knn_solution;

print gmm knn_solution;

diffs = abs(gmm[,2] - knn_solution[,2]);
diff_locs = loc(diffs > 0)`;
print diff_locs;
diff_data = knn[diff_locs, 1];
diff_densities = density[diff_locs];
print diff_data diff_densities;

data = knn[,1];
knn_sol = knn[,2];
gmm_sol = knn[,2];

create plotdata var {diff_data, diff_densities, density, data, knn_sol, gmm_sol};
append;

quit;


/*SORTING THE DATA FOR PLOTS*/
proc sort data = plotdata out = data_sorted;
	by data;
run;


/*The modes plotted on the density*/
proc sgplot data = data_sorted noautolegend;
	series x = data y = density / name = "f(x)"; 
	scatter x = diff_data y = diff_densities /name = "Local Modes" markerattrs =  (symbol = circlefilled color = red);
	title "Differences in cluster assignments between EM and kNN";
	title2 "(observations which received different classifications are marked in red)";
	xaxis label = "x";
	yaxis label = "Density Estimates";
run;
