alpha = csvread('alpha.csv',1,0);
Th = csvread('Th.csv',1,0);
var = csvread('var.csv',1,0);
integ(alpha,Th,var);
csvwrite('ans.csv',ans);
