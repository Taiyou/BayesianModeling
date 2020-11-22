# 波平問題:
# 過去10回の釣果 c(0,1,0,0,2,0,1,0,0,1)
# 事前分布はガンマ分布f(θ|α=6, λ=3)
# 尤度はポアソン分布に従っているとする.
# 奥さんは主人のθ=2, std=0.8
# 
# 6) 事後分布であるガンマ分布 f(θ|α=11, λ=13)を独立MH法を
# シミュレート
# 理論値: 0.846
#
# 7) 事後分布がガウス分布 平均3.0, 分散0.5を独立MH法を用いてシミュレート
# 理論値: 0.846
#
# 8) 波平釣果問題の事後分布であるガンマ分布をランダムウォークMH法でシミュレートしています。10万回サンプリングし, 初めの1000個のバーンイン期間と見積もります。分散を1.0と0.001のふた通りで、トレースラインを観察し、平均を経緯さんして理論値と比較しなさい.
#
#
# by Hiro Taiyo Hamada, Araya Inc.

set.seed(1234);

# initial parameters
n_sample = c(10,100,1000,10000,100000);
burn = 1000;

data = c(0,1,0,0,2,0,1,0,0,1);

# 1) proposal distribution plot
x = -100:500;
x = x/100;
alpha = 11; lambda = 13;
mu = 3; var = 0.5
split.screen(c(2,1)); 
proposal_gammadist = dgamma(x,shape=alpha,rate=lambda);
screen(1); plot(proposal_gammadist, main="Gamma Dist(α=11, λ=13)");
proposal_gaussdist = dnorm(x, mean=mu,sd=sqrt(var));
screen(2); plot(proposal_gaussdist, main="Gauss Dist(μ=3, std=sqrt(0.5))");

# 2) ガンマ分布のカーネル
# α=11, λ=13
gamma_kernel <- function(theta,data){
	if (theta < 0) {
		return(0)
	}
	likelihood = prod(dpois(data, lambda=theta))
	prior      = dgamma(theta, shape=6,rate=3);
	return(prior * likelihood)
}

# independent_MH
independent_MH <- function(data,mu_proposal=1,scale_proposal=sqrt(0.5),theta_init=2,nrun=10) {

	# initial values
	theta_current = theta_init;
	posterior = c(theta_current)
	count = 0; # count accepted proposals

	for (ii in 1:nrun){
		theta_proposal = rnorm(1, mean=mu_proposal, sd=scale_proposal);

		# likelihood x prior
		p_current = gamma_kernel(theta_current,data)
		p_proposal= gamma_kernel(theta_proposal,data)

		# probability to accept proposal
		p_accept =  dnorm(x=theta_current,mean=mu_proposal, sd=scale_proposal) * p_proposal /
        (dnorm(x=theta_proposal,mean=mu_proposal, sd=scale_proposal) * p_current)

		# probability to accept proposal
		if (runif(1) < p_accept){
			count = count + 1;
			theta_current = theta_proposal
		}
		# append posterior
		posterior = append(posterior, theta_current);
	}
	print("Accpeted rate:")
	print(count/nrun*100);
	print("%")
	return(posterior)
}

# 5)問題
posterior = independent_MH(data, nrun=10);
print(mean(posterior))
posterior = independent_MH(data, nrun=100);
print(mean(posterior))
posterior = independent_MH(data, nrun=1000);
print(mean(posterior))
posterior = independent_MH(data, nrun=10000);
print(mean(posterior))

# 6)問題
posterior = independent_MH(data, mu_proposal=3,scale_proposal=sqrt(0.5), nrun=100000);
print(mean(posterior[10000:100000]))


# random_MH
independent_MH <- function(data,scale_proposal=sqrt(0.5),theta_init=4,nrun=10) {

	# initial values
	theta_current = theta_init;
	posterior = c(theta_current)
	count = 0; # count accepted proposals

	for (ii in 1:nrun){
		theta_proposal = rnorm(1, mean=theta_current, sd=scale_proposal);

		# likelihood x prior
		p_current = gamma_kernel(theta_current,data)
		p_proposal= gamma_kernel(theta_proposal,data)

		# probability to accept proposal
		p_accept =  p_proposal/p_current

		# probability to accept proposal
		if (runif(1) < p_accept){
			count = count + 1;
			theta_current = theta_proposal
		}
		# append posterior
		posterior = append(posterior, theta_current);
	}
	print("Accpeted rate:")
	print(count/nrun*100);
	print("%")
	return(posterior)
}

# 8 波平釣果問題
# var = 1.0, 受容率=29.55%, 0.8469911
split.screen(c(2,1)); 
posterior = independent_MH(data, scale_proposal=1,nrun=100000);
print(mean(posterior[10000:100000]));
screen(1);
plot(posterior[1:1000],main="var=1.0");

# var = 0.001, 受容率=95.35%, 0.8375484
posterior = independent_MH(data, scale_proposal=sqrt(0.001),nrun=100000);
print(mean(posterior[10000:100000]));
screen(2); 
plot(posterior[1:1000],main="var=0.001");

