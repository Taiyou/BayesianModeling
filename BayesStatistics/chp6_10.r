# 10) p=10.2, q=5.8のベータ分布をHMC法によってシミュレートし, EAPと理論値である0.638を比較しなさい
# T = 1000, L = 100, ep = 0.01, theta[1] = 0.5
#
#
#
# by Hiro Taiyo Hamada, Araya Inc.

# Leap Frog Method
LFM <- function(theta_prev, p_prev, p_=10.2,q_=5.8,ep=0.01, L=100){
	# leap frog method
	for (ii in 1:L){
		h_theta1   = (1-p_)/theta_prev + (q_-1)/(1-theta_prev)
		p_ep2      = p_prev - ep/2*(h_theta1);
		theta_pred = theta_prev + ep*p_ep2;
		h_theta2   = (1-p_)/theta_pred + (q_-1)/(1-theta_pred)
		p_pred     = p_ep2 - ep/2*h_theta2;

		theta_prev = theta_pred;
		p_prev     = p_pred;
	}
	return(list(p_pred, theta_pred))
}

# potential:
H <- function(theta, p, p_=10.2, q_=5.8){
	h = -(p_-1)*log(theta) - (q_-1)*log(1-theta);
	potential = h + 1/2*p^2;;
	return(potential)
}

# Hamilton Mente Calro
HMC <- function(theta_init=0.5, T=1000){
	# 初期値をおく
	theta    = matrix(0, nrow=T, ncol=1);
	theta[1] = theta_init;
	count    = 0;
	for (ii in 1:T){
		# 独立な標準正規乱数 p(t)を発生させる
		p_prev     = rnorm(1, mean=0, sd=1);

		# LFM法で遷移させ候補点を入手.
		theta_prev = theta[ii];
		prediction = LFM(theta_prev, p_prev);
		p_pred     = prediction[[1]];
		theta_pred = prediction[[2]];

		# 補正係数
		r_coeff = exp(H(theta_prev, p_prev)) - H(theta_pred, p_pred);

		# 受容判定
		if (runif(1) < r_coeff){
			theta[ii+1]   = theta_pred;
			count         = count + 1;
		}else theta[ii+1] = theta_prev;
	}

	acceptance_rate = count/T*100
	return(list(theta,acceptance_rate))
}

results = HMC();
theta           = results[[1]]
acceptance_rate = results[[2]]
sprintf("受容率: %5.1f ", acceptance_rate);
sprintf("推定値: %5.4f ", mean(theta[100:1000]));
sprintf("理論値: %5.4f ", 0.638);
