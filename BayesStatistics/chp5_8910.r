# 基礎からの統計学
# 5章 章末問題
# 8)「正選手問題」の事後分布はp=10.2, q=5.8のベータ分布でした。
# 初期位置は、theta(1)=0.1とし、最初は静止しているものとして,p(1)=0とします
#定数をepisilon=0.05, L=5としてリープフロッグ法として実行し挙動を示せ.
#
#
# 9) episilon=0.01, L=15の場合でも示せ.
#
#
# by Hiro Taiyo Hamada, Araya Inc.

# initial parameters
theta_init = 0.1;
p_init = 0;
epsilon = 0.05;
L = 5;

% Leap-frog method k,
% t, p, theta, h(theta), H(theta, p)を表示する.
% 表5.1 p.116ページ
LFM_51 <- function(theta_init=0.1,p_init=0,epsilon=0.05,L=15){

	# t, p, theta, h(theta), H(theta, p)を保存する.
	alpha    = 11;
	lambda   = 13;
	p        = matrix(0, nrow=L+1, ncol=1);
	theta    = matrix(0, nrow=L+1, ncol=1);
	h        = matrix(0, nrow=L+1, ncol=1);
	H        = matrix(0, nrow=L+1, ncol=1);
	results  = matrix(0, nrow=L, ncol=5);

	# initial values
	p[1]     = p_init;
	theta[1] = theta_init;

	#計算するleap frogを計算する
	for (ii in 1:L){
		# count time
		results[ii, 1] = ii;

		h[ii] = lambda*theta[ii] - (alpha-1)*log(theta[ii]);
		H[ii] = h[ii] + 1/2*p[ii]^2;

		# leap frog method
		p_ep2       = p[ii] - epsilon/2*(lambda-(alpha-1)/theta[ii]);
		theta[ii+1] = theta[ii] + epsilon*p_ep2;
		p[ii+1]     = p_ep2 - epsilon/2*(lambda-(alpha-1)/theta[ii+1])

	}

	results[,2] = p[1:L];
	results[,3] = theta[1:L];
	results[,4] = h[1:L];
	results[,5] = H[1:L];
	return(results)
}

results = LFM()
print(results)


# 章末問題8) 9)
LFM_89 <- function(theta_init=0.1,p_init=0,epsilon=0.05,L=5){

	# t, p, theta, h(theta), H(theta, p)を保存する.
	p_       = 10.2;
	q_       = 5.8;
	p        = matrix(0, nrow=L+1, ncol=1);
	theta    = matrix(0, nrow=L+1, ncol=1);
	h        = matrix(0, nrow=L+1, ncol=1);
	H        = matrix(0, nrow=L+1, ncol=1);
	results  = matrix(0, nrow=L, ncol=5);

	# initial values
	p[1]     = p_init;
	theta[1] = theta_init;

	#計算するleap frogを計算する
	for (ii in 1:L){
		# count time
		results[ii, 1] = ii;

		h[ii] = -(p_-1)*log(theta[ii]) - (q_-1)*log(1-theta[ii]);
		H[ii] = h[ii] + 1/2*p[ii]^2;

		# leap frog method
		h_theta1    = (1-p_)/theta[ii] + (q_-1)/(1-theta[ii])
		p_ep2       = p[ii] - epsilon/2*(h_theta1);
		theta[ii+1] = theta[ii] + epsilon*p_ep2;
		h_theta2    = (1-p_)/theta[ii+1] + (q_-1)/(1-theta[ii+1])
		p[ii+1]     = p_ep2 - epsilon/2*h_theta2;

	}

	results[,2] = p[1:L];
	results[,3] = theta[1:L];
	results[,4] = h[1:L];
	results[,5] = H[1:L];
	print(results)
	return(results)
}

# 章末問題8)
results = LFM_89();

# 章末問題9)
results = LFM_89(theta_init=0.1,p_init=0,epsilon=0.01,L=15)