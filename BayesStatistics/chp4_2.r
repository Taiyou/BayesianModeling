# This is the answer for 2) in chapter 4 problems
#
#
#
# by Hiro Taiyo Hamada, Araya Inc.

n <- c(10, 100, 1000, 10000, 100000, 1000000);
results = 1:length(n);

for (j in 1:length(n)){

	x <- numeric(n[j]);

	for (i in 1:n[j]){ 
		x[i] <- (runif(1)^2+runif(1)^2 <1 )*4;
	};

	results[j] = mean(x)
	print(mean(x));
}