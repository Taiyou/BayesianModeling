# This is the answer for 4) in chapter 4 problems
#
#
#
#
# by Hiro Taiyo Hamada, Araya Inc.

p = matrix(c(0.3, 0.2, 0.5), 1,3);
z <- rbind(c(0.2, 0.2, 0.6), c(0.1, 0.6, 0.3), c(0.3, 0.5, 0.2));
results = matrix(0, nrow=10, ncol=3);

for (i in 1:10){
	p = p %*% x
	results[i,] = p %*% x;
}

print(results)