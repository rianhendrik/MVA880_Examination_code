proc iml;
use sasuser.pcomp; read all into d;
use sasuser.dat2; read all into d2;
use sasuser.dat3; read all into d3;

/*EXPLANATION OF MODULE INPUTS*/
/*dmat - data matrix*/
/*plots - TRUE/FALSE (depending on if you want plots generated or not*/
/*corr - TRUE/FALSE (depending on if you want eigenvectors and values derived from correlation or var-cov matrix*/

start PCA; *(dmat, plots, corr, print);
mean = dmat[:,]; dmatc = dmat-mean; *centering your data; *Should I not maybe standadise it?;


n = nrow(dmat); p = ncol(dmat); *defining row and column dimensions of the data set;
covmat = (dmatc`*dmatc)/(n-1); *calculating variance-covariance matrix of data;
stdcol = sqrt(dmatc[##,]/(n-1)); stddmatc = dmatc/stdcol; corrmat = (stddmatc`*stddmatc)/(n-1); *calculating correlation matrix of data;

if corr = "true" then call Eigen(evalues, evectors, corrmat); *computing eigenvalues and vectors using corrmat;
else call Eigen(evalues, evectors, covmat); *computing eigenvalues and vectors using covmat;
* print evectors,,evalues;
pcomps = dmatc*evectors; *computing p principle components;

/* graph the simulated data */
title "Determining the right number of principle components to choose";
/*Skree plot*/
ods graphics reset = all;
x = do(1, p, 1)`; y = evalues;
create Skree var {x y}; append; close Skree; /* write data */
if plots = "true" then submit;
  proc sgplot data = Skree;
  	title "Skree plot";
    series x = x y = y;
   	xaxis grid label="Principal component";
    yaxis grid label="Eigenvalues";
  run;
endsubmit;

/*Variance plot*/
den = sum(evalues); *denominator of variance explained value;
varex_cum = J(p, 1, .); varex = J(p, 1, .);
do i = 1 to p;
	varex_cum[i] = sum(evalues[1:i])/den; *prop variance explained by 1:ith component;
	varex[i] = evalues[i]/den;
end;

create varexplot var {x varex varex_cum}; append; close varexplot;
if plots = "true" then submit;
	proc sgplot data = varexplot;
		title "Proportion of variance explained by Principle components (individual and cumulative)";
		series x = x y = varex_cum / lineattrs=(pattern=mediumdash color=blue thickness=1);
		series x = x y = varex;
		xaxis grid label="Principle component";
		yaxis grid label="Proportion of variance explained";
	run;
endsubmit;

title "Determining the right number of principle components to choose";
g1 = ncol(loc(evalues>1));
print "The number of components with eigenvalues greater then 1 is" g1;

/*Determining the correlation between each component and each original data variable*/

evalmat = repeat(evalues, ncol(evalues), nrow(evalues));
row_titles = do(1, p, 1);
col_titles = do(0, p, 1);
influence = col_titles//(t(row_titles)||(evectors#sqrtevalmat)));

finish;

dmat=d2; corr="true"; plots="true";

call PCA; *(d2, "true", "true", "true");

print evalues evectors,, influence;




