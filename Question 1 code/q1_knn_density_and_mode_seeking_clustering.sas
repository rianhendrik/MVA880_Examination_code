ods html close;
ods html;

proc iml;
use sasuser.exam_dist;
read all into x;

/*THE kNN SUBROUTINE*/

start kNN_mode_clus;        


/*ESTIMATING THE DENSITY*/
n = nrow(data);
p = ncol(data);
*serials = (1:n)`; *since we are sampling, these serials must not reset the actual serials attached to data points in the sample.;
knn_serials = (1:n)`;
d = distance(data);
fx_h = J(n, 1, .);
do e = 1 to n;
	d_i = d[, e];
	call sort(d_i);
	fx_h[e] = 1/d_i[k];
end;

control_data = data||knn_serials||fx_h;



/*GRADIENT ASCENT*/
modes_knn = J(n, 1, .); /*Each x value will have a corresponding serial (index) number of the mode to which it belongs*/
modes_true = J(n, 1, .);
fx_h_modes = J(n, 1, .); /*We also want the density values at these modes, to plot them over the density estimate*/
do e = 1 to n;
	mat = d[, e]||fx_h||knn_serials||data[, p];
	call sort(mat, {1});
	maxi = 0;
	window = mat[1:k, ];
	loc_max_f = window[<:>, 2];
	mode_x_knn = mat[loc_max_f, 3];
	true_mode_x = mat[loc_max_f, 4];
	do while (mode_x_knn ^= maxi);
		maxi = mode_x_knn; *This is the serial of the local max x value;
		mat2 = d[, maxi]||fx_h||knn_serials||data[, p];
		call sort(mat2, {1});
		window2 = mat2[1:k, ];
		loc_max_f2 = window2[<:>, 2];
		mode_x_knn = mat2[loc_max_f2, 3];
		true_mode_x = mat2[loc_max_f2, 4];
	end;
	modes_knn[e] = mode_x_knn; *modes according to the knn_serials;
	modes_true[e] = true_mode_x;
	fx_h_modes[e] = mat2[loc_max_f2, 2];
end;
results = data||modes_true;

finish;


data = x;
k = 30;
call kNN_mode_clus;

knn = results;

modes = unique(knn[,2]);
print modes;
clusts = J(nrow(knn), 1, .);
do i = 1 to nrow(knn);
	do j = 1 to 3;
		if knn[i, 2] = modes[j] then clusts[i] = j;
	end;
end;

knn_solution = knn[,1]||clusts; print knn_solution;

mus_knn = modes`;
stds_knn = J(3, 1, .);
pis_knn = J(3, 1, .);
do i = 1 to 3;
	clust_i = knn_solution[loc(knn_solution[,2] = i),1];
	stds_knn[i] = std(clust_i);
	pis_knn[i] = nrow(clust_i)/nrow(knn_solution);
end;

print mus_knn stds_knn pis_knn;
