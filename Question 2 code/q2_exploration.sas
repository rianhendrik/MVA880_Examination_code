ods html close;
ods html;

data true_digits;
set sasuser.exam_digits;
keep y;
run;


proc freq data = true_digits;
run;



proc iml;
use sasuser.exam_digits;
read all into digits;


ods listing gpath = "C:\Users\rianh\Dropbox\MVA880\EXAM\LaTeX document\Exam document\figures"
image_dpi=300;


true_d = digits[,1];
d = digits[,2:ncol(digits)];

ndigits = nrow(d);
npixels = ncol(d);
head = d[1:5,]; print head;
print ndigits npixels;

i = 20;
digit_i = shape(d[i,], 16, 16);
call HeatmapCont(digit_i);
i = 30;
digit_i = shape(d[i,], 16, 16);
call HeatmapCont(digit_i);


quit;

