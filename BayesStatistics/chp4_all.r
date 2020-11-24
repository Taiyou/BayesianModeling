#
#
#
#
#
# by Hiro Taiyo Hamada, Araya Inc.


library(ggplot2)

d <- read.csv(file='input/data-salary.txt')
p <- ggplot(data=d, aes(x=X, y=Y))
p <- p + theme_bw(base_size=18)
p <- p + geom_point(shape=1, size=3)
ggsave(file='output/fig4-2.png', plot=p, dpi=300, w=4, h=3)

# 年齢が一つ上がるごとに21.9万円上がることになる.
res_lm <- lm(Y ~ X, data=d)
res_lm

# 4.4.4 Rのlm関数で推定
# lm関数の結果を使うと,パラメータの信頼区間と新しいデータの
# 予測区間を求めることができる.
X_new <- data.frame(X=23:60)
conf_95 <- predict(res_lm, X_new, interval='confidence', level=0.95)
conf_95 <- data.frame(X_new, conf_95)
conf_50 <- predict(res_lm, X_new, interval='confidence', level=0.50)
conf_50 <- data.frame(X_new, conf_50)

pred_95 <- predict(res_lm, X_new, interval='prediction', level=0.95)
pred_95 <- data.frame(X_new, pred_95)
pred_50 <- predict(res_lm, X_new, interval='prediction', level=0.50)
pred_50 <- data.frame(X_new, pred_50)

# visualize confidence
p <- ggplot()
p <- p + theme_bw(base_size=18)
# plotting confidence interval
p <- p + geom_ribbon(data=conf_50, aes(x=X, ymin=lwr, ymax=upr), alpha=1/6)
p <- p + geom_ribbon(data=conf_95, aes(x=X, ymin=lwr, ymax=upr), alpha=2/6)
# plotting line
p <- p + geom_line(data=conf_50, aes(x=X, y=fit), size=1)
# plotting datapoints
p <- p + geom_point(data=d, aes(x=X, y=Y), shape=1, size=3)
# plotting coordinates
p <- p + labs(x='X', y='Y') + coord_cartesian(xlim=c(22, 61), ylim=c(200, 1400))
p <- p + scale_y_continuous(breaks=seq(from=200, to=1400, by=400))
p
# save figure.
# ggsave(p, file='output/fig4-3-left.png', dpi=300, w=4, h=3)

# visualize prediction
p <- ggplot()
p <- p + theme_bw(base_size=18)
# plotting confidence interval
p <- p + geom_ribbon(data=pred_50, aes(x=X, ymin=lwr, ymax=upr), alpha=1/6)
p <- p + geom_ribbon(data=pred_95, aes(x=X, ymin=lwr, ymax=upr), alpha=2/6)
# plotting line
p <- p + geom_line(data=pred_50, aes(x=X, y=fit), size=1)
# plotting datapoints
p <- p + geom_point(data=d, aes(x=X, y=Y), shape=1, size=3)
# plotting coordinates
p <- p + labs(x='X', y='Y') + coord_cartesian(xlim=c(22, 61), ylim=c(200, 1400))
p <- p + scale_y_continuous(breaks=seq(from=200, to=1400, by=400))
p
# save figure.
# ggsave(p, file='output/fig4-3-left.png', dpi=300, w=4, h=3)


# 4.4.5 Stanで実装
setwd("Desktop/Bayesstatistical/RStanBook-master/chap04")
library(rstan)

d <- read.csv(file='input/data-salary.txt')
data <- list(N=nrow(d), X=d$X, Y=d$Y)
fit  <- stan(file='model4-5.stan', data=data, seed=1234)

# 4.4.8 収束診断をファイルへ出力する
# MCMCの収束判定パッケージ
install.packages("ggmcmc")

# rstan-save-diagnostics.r
library(rstan)
library(ggmcmc)

load('output/result-model4-5.Rdata')

write.table(data.frame(summary(fit)$summary), 
	file='output/fit-summary.txt', sep='\t', quote=FALSE, col.name=NA)

ggmcmc(ggc(fit, inc_warmup=TRUE, stan_include_auxiliar=TRUE),
	file='output/fit-traceplot.pdf', plot='traceplot')

ggmcmc(ggs(fit), file='output/fit-ggmcmc.pdf')


# 4.4.9 MCMCの設定の変更
# warmup期間を短くするなど, MCMCの設定を変えたいことは頻繁にある.
library(rstan)

d <- read.csv(file='input/data-salary.txt')
data <- list(N = nrow(d), X=d$X, Y=d$Y)


# 4.4.10
# 並列計算の実行方法
# stan関数などの前に以下のRコードを書き込む
# rstan_options(auto_write=TRUE)
# options(mc.cores=parallel::detectCores())

stanmodel <- stan_model(file='model/model4-5.stan')

fit <- sampling(
	stanmodel,
	data = data,
	pars = c('b', 'sigma')
	# 初期値の設定
	init = function(){
		list(a=runif(1,-10,10), b=runif(1,0,10), sigma=10)
	},
	seed = 123,
	# iterationのステップ数などを入れている
	# chanis: 最低3は欲しい
	# iter: 最終的には1000~1500
	# warmupはtrace plotを使ってみることが多いが, 100~500が多い
	# thin: 通常は1で実行している.
	chains = 3, iter = 1000, warmup = 200, thin = 2
	)

# 4.1.11
# ベイズ信頼区間とベイズ予測区間の算出
# rstan-extract-MCMsample.r
load('output/result-model4-5.RData')

ms <- rstan::extract(fit)

# 傾きbのMCMCサンプルをR Console上で確認し,
ms$b
# 95%ベイズ信頼区間を求めてみよう.
quantile(ms$b, probs=c(0.025, 0.975))

# 4.4.12


head(ms$y_new[ ,1:4])


