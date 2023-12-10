## -----------------------------------------------------------------------------
IPW <- function(M, K=2) {
  n <- nrow(M)
  d <- ncol(M)
  Y <- ifelse(M == 0, 0, 1)  # 使用向量化操作代替循环
  
  p_hat <- sum(Y) / (n * d)
  
  Sigma_phat <- crossprod(M) - (1 - p_hat) * diag(colSums(M^2))
  Sigma_tphat <- crossprod(t(M)) - (1 - p_hat) * diag(rowSums(M^2))

  
  eigen_Sigma_phat <- eigen(Sigma_phat)
  eigen_Sigma_tphat <- eigen(Sigma_tphat)
  
  # 得到分解
  V <- eigen_Sigma_phat$vectors[, 1:K]
  lambda_hat <- sqrt(eigen_Sigma_phat$values[1:K] - sum(eigen_Sigma_phat$values[(K+1):d]) / (d - K -1)) / p_hat
  U <- eigen_Sigma_tphat$vectors[, 1:K]
  
  return(list(V = V, U = U, lambda_hat = lambda_hat))
}

## -----------------------------------------------------------------------------
initial<-function(n,d){
  A <- matrix(runif(n * 2, min = -5, max = 5), n, 2)
  B <- matrix(runif(d * 2, min = -5, max = 5), d, 2)
  M0 <- A%*%t(B)
  
  return(M0)
}

## -----------------------------------------------------------------------------
#xtable::xtable(IPW(initial(n = 10, d = 2))$V)

## ----eval = FALSE-------------------------------------------------------------
#  IPW(initial(n = 10, d = 2))$U

## -----------------------------------------------------------------------------
my_sample <- function(values, n, probabilities = rep(1/length(values), length(values))) {
  # Compute the cumulative distribution function
  cdf <- cumsum(probabilities)
  
  # Generate random numbers
  samples <- numeric(n)
  for (i in 1:n) {
    u <- runif(1)
    samples[i] <- values[findInterval(u, cdf)+1]
  }
  
  return(samples)
}

## -----------------------------------------------------------------------------
# Initialize the values to be sampled and the sample size
values <- 1:10
n <- 10

random_samples <- my_sample(values, n)
print(random_samples)

## -----------------------------------------------------------------------------
# Generate random samples from the standard Laplace distribution using the inverse transform method
inverse_sample_laplace <- function(n) {
  u <- runif(n)
  samples <- ifelse(u < 0.5, log(2 * u), -log(2 - 2 * u))
  return(samples)
}

# Generate 1000 random samples from the standard Laplace distribution
n <- 1000
laplace_samples <- inverse_sample_laplace(n)

# Plot the histogram of the generated samples
hist(laplace_samples, freq = FALSE, main = "Samples vs Laplace Distribution",
     xlab = "Samples", ylab = "Density", col = "lightblue")

# Plot the density function curve of the standard Laplace distribution
x <- seq(-5, 5, length.out = 100)
density_laplace <- 0.5 * exp(-abs(x))
lines(x, density_laplace, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Generated Samples", "Laplace Distribution"),
       col = c("lightblue", "red"), lwd = 2, bty = "n")

## -----------------------------------------------------------------------------
# Rejection sampling to generate random samples from the beta distribution
generate_beta_sample <- function(n, a, b) {
  # Let g(x) = 1, 0 < x < 1 and M = max(a, b) + 1. 
  samples <- numeric(n)
  count <- 0
  M <- max(a, b) + 1
  
  while (count < n) {
    u <- runif(1)
    x <- runif(1)
    
    if (dbeta(x, a, b) > M*u) {
      samples[count + 1] <- x
      count <- count + 1
    }
  }
  
  return(samples)
}

## -----------------------------------------------------------------------------
# Generate random samples from the Beta(3, 2) distribution
n <- 1000
a <- 3
b <- 2
beta_samples <- generate_beta_sample(n, a, b)

# Plot the histogram of the generated samples
hist(beta_samples, freq = FALSE, main = "Samples vs. Beta(3,2) Distribution",
     xlab = "Samples", ylab = "Density", col = "lightblue")

# Plot the density function curve of the Beta(3, 2) distribution
x <- seq(0, 1, length.out = 100)
density_beta <- dbeta(x, a, b)
lines(x, density_beta, col = "red", lwd = 2)

# Add a legend
legend("topleft", legend = c("Generated Samples", "Beta(3,2) Distribution"),
       col = c("lightblue", "red"), lwd = 2, bty = "n")

## -----------------------------------------------------------------------------
# Generate samples using the method described in the question
generate_epanechnikov_sample <- function(n) {
  samples <- numeric(n)
  
  for (i in 1:n) {
    U1 <- runif(1, -1, 1)
    U2 <- runif(1, -1, 1)
    U3 <- runif(1, -1, 1)
    
    if (abs(U3) >= abs(U2) && abs(U3) >= abs(U1)) {
      samples[i] <- U2
    } else {
      samples[i] <- U3
    }
  }
  
  return(samples)
}

## -----------------------------------------------------------------------------
# Generate a large number of simulated random samples
n <- 10000
epanechnikov_samples <- generate_epanechnikov_sample(n)

# Construct a histogram density estimate
hist(epanechnikov_samples, prob = TRUE, main = "Histogram Density Estimate",
     xlab = "Samples", ylab = "Density", col = "lightblue")

# Construct a histogram density estimate
x <- seq(-1, 1, length.out = 100)
density_epanechnikov <- 3/4 * (1 - x^2)
lines(x, density_epanechnikov, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Estimate", "Kernel"),
       col = c("lightblue", "red"), lwd = 2, bty = "n")

## -----------------------------------------------------------------------------
# 构造pi_hat的输出函数
var_pi <- function(rho, seed = 1) {
  set.seed(seed)
  # 直接假定d = 1 ，则此时l = rho 
  m <- 1e6
  X <- runif(m, 0, 1 / 2)
  Y <- runif(m, 0, pi / 2)
  pihat <- 2 * rho / mean(rho / 2 * sin(Y) > X)
}

## -----------------------------------------------------------------------------
# 构造函数，使得得到多个pi_hat以计算方差
var_pi_monte <- function(rho, monte) {
  var_pi_result <- sapply(1:monte, function(i) var_pi(rho, i))
  return(round(var(var_pi_result),6))
}

## -----------------------------------------------------------------------------
rho_values <- c(0.5, 0.8, 1)
results <- sapply(rho_values, var_pi_monte, monte = 100)
message(paste("当分别代入", paste(rho_values, collapse = ", "), "时，我们有方差分别为:", paste(results, collapse = ", ")))

## -----------------------------------------------------------------------------
m <- 1e3

# 创建估计函数，输出theta_hat的估计，theta_hat的antithetic估计，和真实的theta
estimate <- function() {
  x <- runif(m, min = 0, max = 1)
  theta_hat <- mean(exp(x))  # Simple Monte Carlo method
  x1 <- runif(m / 2, min = 0, max = 1)
  theta_hat_antithetic <- (mean(exp(x1))+mean(exp(1-x1)))/2  # Antithetic variate approach
  true <- exp(1) - 1
  c(theta_hat, theta_hat_antithetic, true)
}

# 重复以便计算出方差
result <- replicate(1000, estimate())
var_simple <- var(result[1, ])
var_antithetic <- var(result[2, ])
percent_reduction_var <- (var_simple - var_antithetic) / var_simple

# 输出方差及其减少百分比
output <- paste(" Simple Monte Carlo method:", var_simple, "\n",
                "Antithetic variate approach:", var_antithetic, "\n",
                "Empirical estimate of the percent reduction:", percent_reduction_var, "\n",
                "True percent reduction:", (exp(2) - 2 * exp(1) - 1) / (-2 * exp(2) + 8 * exp(1) - 6))

cat(output)

## -----------------------------------------------------------------------------
data <- c(2.14, 2.10, 2.13, 2.15, 2.13, 2.12, 2.13, 2.10, 
          2.15, 2.12, 2.14, 2.10, 2.13, 2.11, 2.14, 2.11)

var(data)

## ----message=FALSE------------------------------------------------------------
library(bayesmeta)  # 引入 bayesmeta 包，用于生成瑞利分布

## -----------------------------------------------------------------------------

# 定义函数 g(x)
g <- function(x) {
  x^2 / sqrt(2*pi) * exp(-x^2/2)
}

xs <- seq(0, 10, 0.1)  # 生成从 0 到 10，步长为 0.1 的序列作为 x 值

# 计算函数 g(x) 在 xs 上的取值
ys.g <- g(xs)

# 生成瑞利分布在 xs 上的取值，其中 scale 参数为 1.5
ys.rayleigh <- drayleigh(xs, scale = 1.5)

# 生成正态分布在 xs 上的取值，其中 mean 参数为 1.5
ys.norm <- dnorm(xs, mean = 1.5)

lim <- max(c(ys.g, ys.rayleigh, ys.norm))  # 计算三个函数的取值的最大值

# 绘制函数曲线
plot(xs, ys.g, type = "l", ylim = c(0, lim))  # 绘制函数 g(x) 的曲线，y 轴范围为 0 到 lim
lines(xs, ys.rayleigh, col = "red", ylim = c(0, lim))  # 绘制瑞利分布的曲线，颜色为红色，y 轴范围为 0 到 lim
lines(xs, ys.norm, col = "blue", ylim = c(0, lim))  # 绘制正态分布的曲线，颜色为蓝色，y 轴范围为 0 到 lim

## -----------------------------------------------------------------------------
M <- 1e4  # 总样本量

# 定义函数 g(x)
g <- function(x) {
  x^2 / sqrt(2*pi) * exp(-x^2/2) * (x > 1)
}

# 定义函数 f1(x)
f1 <- function(x) {
  drayleigh(x, scale = 1.5) * (x > 1)
}

# 定义函数 f2(x)
f2 <- function(x) {
  dnorm(x, mean = 1.5) * (x > 1)
}

# 定义函数 rf1()，生成服从瑞利分布的随机样本
rf1 <- function() {
  rrayleigh(M, scale = 1.5)
}

# 定义函数 rf2()，生成服从正态分布的随机样本
rf2 <- function() {
  rnorm(M, mean = 1.5)
}

# 定义函数 is.rayleigh()，估计瑞利分布的参数 theta1
is.rayleigh <- function() {
  xs <- rf1()  # 生成瑞利分布的随机样本
  mean(g(xs)/f1(xs), na.rm = TRUE)  # 计算估计值，取平均并忽略缺失值
}

# 定义函数 is.norm()，估计正态分布的参数 theta2
is.norm <- function() {
  xs <- rf2()  # 生成正态分布的随机样本
  mean(g(xs)/f2(xs), na.rm = TRUE)  # 计算估计值，取平均并忽略缺失值
}

theta1 <- is.rayleigh()  # 估计瑞利分布的参数 theta1
theta2 <- is.norm()  # 估计正态分布的参数 theta2

theta1  # 打印估计的瑞利分布参数 theta1
theta2  # 打印估计的正态分布参数 theta2

## -----------------------------------------------------------------------------
M <- 1e4  # 总样本量
N <- 1000  # 重复次数
k <- 5  # 分层数

# 逆变换方法生成密度函数为 f_k(x) 的随机数
# 注意 (a, b) 对应不同层的密度函数
inv_fun <- function(n, a, b) {  
  u <- runif(n)
  x <- -log(exp(-a) - (exp(-a) - exp(-b)) * u)
  return(x)
} 

# 使用 sapply 函数进行重复计算
res3 <- sapply(1:N, FUN = function(o){
  x <- inv_fun(M, 0, 1)
  M1 <- mean((1 - exp(-1)) / (1 + x^2))  # 重要抽样均值
  M2 <- numeric(k)
  for (j in 0:(k-1)){
    a <- j / k
    b <- (j + 1) / k
    xj <- inv_fun(M / k, a, b)  # 生成分层随机数
    M2[j+1] <- mean((exp(-a) - exp(-b)) / (1 + xj^2))  # 分层重要抽样均值
  }
  c(M1, sum(M2))
})

# 计算估计方差
c(var(res3[1,]), var(res3[2,]))  # 使用 SSE 方法估计方差

## -----------------------------------------------------------------------------
exercise_6_4 <- function(seed = 123) {
  set.seed(seed)
  n <- 20  # 样本量
  alpha <- 0.05  # 置信水平
  m <- 1000  # Monte Carlo 实验次数

  cv.t <- sapply(1:m, FUN = function(o) {
    y <- rnorm(n)
    c <- qt(1 - alpha/2, n - 1)  # t 分布的上 α/2 分位数
    m <- mean(y)  # 均值估计
    se <- sqrt(var(y))  # 标准误差估计
    as.numeric((m - c * se / sqrt(n) < 0) & (m + c * se / sqrt(n) > 0))  # 置信区间
  })

  level <- mean(cv.t)  # Monte Carlo 实验的均值

  return(data.frame(level = level))
}

exercise_6_4()

## -----------------------------------------------------------------------------
exercise_6_5 <- function(seed = 123) {
  set.seed(seed)
  n <- 20
  t_value <- qt(0.975, n - 1)  # 0.975分位数的t-分布值，用于计算置信区间
  m <- 1000
  cv.t <- sapply(1:m, FUN = function(o) {
    x <- rchisq(n, 2)  # 从自由度为2的卡方分布中抽取样本（注意这里的x是从卡方分布取样）
    mean_x <- mean(x)  # 样本均值的估计值
    se <- sqrt(var(x))  # 样本标准误差的估计值
    as.numeric((mean_x - t_value * se / sqrt(n) < 2) & (mean_x + t_value * se / sqrt(n) > 2))  # 计算置信区间
  })
  level1 <- mean(cv.t)  # Monte Carlo实验的均值，用于估计概率

  # 我们可以得出概率不等于0.95，小于0.95。在example6.4中使用卡方分布来估计方差（真值为4）
  alpha <- 0.05
  UCL <- replicate(1000, expr = {
    x <- rchisq(n, 2)
    (n - 1) * var(x) / qchisq(alpha, df = n - 1)
  })
  # 计算包含sigma^2=4的区间数
  level2 <- sum(UCL > 4) / m
  
  return(data.frame(level1, level2))
}

exercise_6_5(1012)

## -----------------------------------------------------------------------------
exercise_6_A <- function(seed) {
  set.seed(seed)
  num <- c(50, 100, 200, 500, 1000)  # 不同样本量
  m <- 10000  # Monte Carlo 实验次数
 
  er <- NULL
  for (n in num) {
    cv <- qt(0.975, n - 1)  # t 分布的上 0.975 分位数
    
    # 估计卡方分布的第一类错误
    er1 <- mean(sapply(1:m, FUN = function(o) {
      x <- rchisq(n, 1)
      m <- mean(x)
      se <- sqrt(var(x))
      abs((m - 1) * sqrt(n) / se) >= cv
    }))
    
    # 估计均匀分布的第一类错误
    er2 <- mean(sapply(1:m, FUN = function(o) {
      x <- runif(n, 0, 2)
      m <- mean(x)
      se <- sqrt(var(x))
      abs((m - 1) * sqrt(n) / se) >= cv
    }))
    
    # 估计指数分布的第一类错误
    er3 <- mean(sapply(1:m, FUN = function(o) {
      x <- rexp(n, 1)
      m <- mean(x)
      se <- sqrt(var(x))
      abs((m - 1) * sqrt(n) / se) >= cv 
    }))
    
    er <- cbind(er, c(er1, er2, er3))
  }
  
  colnames(er) <- num
  rownames(er) <- c("chi(1)", "U(0,2)", "exp(1)")
  
  return(er)                
}

exercise_6_A(1012)

## ----message=FALSE------------------------------------------------------------
library(purrr)
library(bootstrap)
library(boot)

## -----------------------------------------------------------------------------
# 设置随机种子以确保结果可重复
set.seed(123)

# 定义参数
m <- 1000
alpha <- 0.1

# 进行模拟
simulations <- replicate(m, {
  # 生成m个p值
  p_values <- c(runif(m * 0.95), rbeta(m * 0.05, 0.1, 1))
  
  # Bonferroni校正
  bonferroni_p <- p.adjust(p_values, method = "bonferroni")
  
  # Benjamini-Hochberg校正
  bh_p <- p.adjust(p_values, method = "BH")
  
  # 根据阈值确定是否拒绝原假设
  bonferroni_reject <- bonferroni_p <= alpha
  bh_reject <- bh_p <= alpha
  
  # 计算FWER、FDR和TPR
  bonferroni_fwer <- mean(bonferroni_reject)
  bonferroni_fdr <- mean(bonferroni_reject) / sum(bonferroni_reject)
  bonferroni_tpr <- sum(bonferroni_reject[1:(m * 0.95)]) / (m * 0.95)
  
  bh_fwer <- mean(bh_reject)
  bh_fdr <- mean(bh_reject) / sum(bh_reject)
  bh_tpr <- sum(bh_reject[1:(m * 0.95)]) / (m * 0.95)
  
  # 返回结果
  c(bonferroni_fwer, bonferroni_fdr, bonferroni_tpr, bh_fwer, bh_fdr, bh_tpr)
})

# 计算平均FWER、FDR和TPR
mean_bonferroni_fwer <- mean(simulations[1, ])
mean_bonferroni_fdr <- mean(simulations[2, ])
mean_bonferroni_tpr <- mean(simulations[3, ])
mean_bh_fwer <- mean(simulations[4, ])
mean_bh_fdr <- mean(simulations[5, ])
mean_bh_tpr <- mean(simulations[6, ])


# 输出结果表格
result_table <- data.frame(
  FWER = c(mean_bonferroni_fwer, mean_bh_fwer),
  FDR = c(mean_bonferroni_fdr, mean_bh_fdr),
  TPR = c( mean_bonferroni_tpr, mean_bh_tpr)
)

rownames(result_table) <- c("Bonferroni", "Benjamini-Hochberg")

print(result_table)

## -----------------------------------------------------------------------------
# 设置λ的真实值
true_lambda <- 2

# 设置样本大小
sample_sizes <- c(5, 10, 20)

# 设置bootstrap重复次数
B <- 1000

# 设置模拟次数
m <- 1000


# 执行bootstrap估计的函数
bootstrap_estimation <- function(sample, B) {
  bootstrap_estimates <- rep(0, B)
  
  for (i in 1:B) {
    bootstrap_sample <- sample(sample, size = length(sample), replace = TRUE)
    bootstrap_estimates[i] <- 1/mean(bootstrap_sample) - 1/mean(sample)
  }
  
  bootstrap_mean <- mean(bootstrap_estimates)
  bootstrap_sd <- sd(bootstrap_estimates)
  
  c(bootstrap_mean = bootstrap_mean, bootstrap_sd = bootstrap_sd)
}

# 对每个样本大小进行模拟研究
results <- lapply(sample_sizes, function(n) {
  replicate(m, {
    # 从指数分布生成样本
    sample <- rexp(n, rate = true_lambda)
    
    # 执行bootstrap估计
    bootstrap_results <- bootstrap_estimation(sample, B)
    
    bootstrap_results
  })
})

# 计算平均bootstrap估计和bootstrap标准差
mean_bootstrap_estimate <- sapply(results, function(x) mean(x[1, ]))
mean_bootstrap_sd <- sapply(results, function(x) mean(x[2, ]))

# 计算理论偏差和标准误差
theoretical_bias <- true_lambda / (sample_sizes - 1)
theoretical_se <- true_lambda * sample_sizes / ((sample_sizes - 1) * sqrt(sample_sizes - 2))

# 创建数据框来存储结果
result_table <- data.frame(
  样本大小 = sample_sizes,
  平均Bootstrap估计 = mean_bootstrap_estimate,
  平均Bootstrap标准差 = mean_bootstrap_sd,
  理论偏差 = theoretical_bias,
  理论标准误差 = theoretical_se
)

print(result_table)

## -----------------------------------------------------------------------------
example_7_3 <- function(seed = 1012) {
  set.seed(seed)
  
  # 定义boot.t.ci函数，用于计算bootstrap t置信区间
  boot.t.ci <- function(x, B = 500, R = 100, level = 0.95, statistic) {
    stat <- numeric(B) # 存储bootstrap样本的统计量
    se <- numeric(B) # 存储bootstrap样本的标准误差
    
    # 定义boot.se函数，用于计算bootstrap样本的标准误差
    boot.se <- function(x, R, f) {
      th <- replicate(R, expr = {
        i <- sample(1:nrow(x), size = nrow(x), replace = TRUE)
        f(x[i, ])
      })
      sd(th) # 返回bootstrap样本的标准差
    }
    
    for (b in 1:B) {
      y <- x[sample(nrow(x), replace = TRUE), ] # 构建bootstrap样本
      stat[b] <- statistic(y) # 计算bootstrap样本的统计量
      se[b] <- boot.se(y, R = R, f = statistic) # 计算bootstrap样本的标准误差
    }
    
    stat0 <- statistic(x) # 原始样本的统计量
    t.stats <- (stat - stat0) / se # 计算t统计量
    se0 <- sd(stat) # 原始样本的标准误差
    alpha <- 1 - level
    Qt <- quantile(t.stats, c(alpha/2, 1 - alpha/2), type = 1) # 计算t统计量的分位数
    CI <- rev(stat0 - Qt * se0) # 构建置信区间
    
    CI
  }
  
  dat <- cbind(law$LSAT, law$GPA) # 原始数据集
  stat <- function(dat) {
    mean(cor(dat[, 1], dat[, 2])) # 统计量函数，计算LSAT和GPA的相关系数的均值
  }
  
  ci <- boot.t.ci(dat, B = 2000, R = 200, statistic = stat) # 执行bootstrap t置信区间计算
  print(ci) # 输出结果
}

example_7_3() # 调用示例函数

## ----message=FALSE------------------------------------------------------------
library(boot) 

## -----------------------------------------------------------------------------
example_7_5 <- function(seed = 1023) {
  set.seed(seed)
  dat <- aircondit$hours  # 获取数据集 aircondit 中的 hours 列作为输入数据
  
  theta.boot <- function(dat, ind) {
    mean(dat[ind])
  }
  
  boot.obj <- boot(dat, statistic = theta.boot, R = 2000)  # 进行自助法计算
  
  print(boot.obj) 
  print(boot.ci(boot.obj, type = c("basic", "norm", "perc", "bca")))  # 输出自助法置信区间
}

example_7_5()


## ----message=FALSE------------------------------------------------------------
  library(bootstrap)

## -----------------------------------------------------------------------------
example_7_8 <- function(seed = 1023) {
  set.seed(seed)

  
  sc <- scor
  
  theta <- function(x) {
    sigma <- cov(x)  # 计算协方差矩阵
    pca.sigma <- prcomp(sigma)  # 进行主成分分析
    theta <- pca.sigma$sdev[1] / sum(pca.sigma$sdev)  # 获取主成分分析的第一个特征值占比
    theta
  }
  
  n <- NROW(sc)
  theta.j <- numeric(n)
  
  # 计算每个样本的 theta 值
  for (i in 1:n) {
    theta.j[i] <- theta(sc[-i, ])
  }
  
  theta.hat <- theta(sc)
  
  bias <- (n - 1) * (mean(theta.j) - theta.hat)  # 计算偏差
  se <- sqrt((n - 1) * var(theta.j))  # 计算标准误差
  
  result <- list(Bias.jack = bias, SE.jack = se)
  return(result)
}

example_7_8()

## ----message=FALSE------------------------------------------------------------
library(DAAG)
attach(ironslag)

## -----------------------------------------------------------------------------
n <- length(magnetic)  # magneti和chemical是数据集中的向量，在DAAG ironslag中
e1 <- e2 <- e3 <- e4 <- numeric(n)

for (k in 1:n) {
  y <- magnetic[-c(k, k+1)]  # 离开两个样本后的目标变量
  x <- chemical[-c(k, k+1)]  # 离开两个样本后的特征变量
  
  J1 <- lm(y ~ x)  # 拟合线性模型
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
  e1[k] <- magnetic[k] - yhat1 
  
  J2 <- lm(y ~ x + I(x^2))  # 拟合二次模型
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
  e2[k] <- magnetic[k] - yhat2 
  
  J3 <- lm(log(y) ~ x)  # 拟合对数线性模型
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] 
  yhat3 <- exp(logyhat3)  
  e3[k] <- magnetic[k] - yhat3 
  
  J4 <- lm(log(y) ~ log(x))  # 拟合对数对数模型
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])  
  yhat4 <- exp(logyhat4) 
  e4[k] <- magnetic[k] - yhat4  
}

# 输出平均平方误差
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## ----message=FALSE------------------------------------------------------------
  library(RVAideMemoire)

## -----------------------------------------------------------------------------
feed1 = "soybean"
feed2 = "linseed"

data <- chickwts

# 提取大豆饲料和亚麻籽饲料的体重数据
xs <- sort(data$weight[data$feed == feed1])[1:12]
ys <- sort(data$weight[data$feed == feed2])[1:12]

rep <- 1000

# 将两组数据合并为一个向量
zs <- c(xs, ys)

n1 <- length(xs)
n2 <- length(ys)
n <- n1 + n2

Ts <- numeric(rep)

# 进行重复次数为rep的循环
for (i in 1:rep) {
  # 从合并的数据中随机选择n1个样本，作为第一组样本
  ks <- sample(1:n, n1, replace = FALSE)
  zs1 <- zs[ks]
  zs2 <- zs[-ks]
  Ts[i] <- cramer.test(zs1, zs2)$statistic
}

# 计算Cramer-von Mises检验的统计量和p值
cvm <- cramer.test(xs, ys)
T.hat <- cvm$statistic
p.hat <- mean(abs(T.hat) < abs(Ts))

hist(Ts)

## -----------------------------------------------------------------------------
n1 = 100
n2 = 200

# 返回H1的概率，其中H0是等方差（即sd1 = sd2）
test_equal_variance = function(xs1, xs2){
  # 让x1小于x2
  if (length(xs2) < length(xs1)) {
    tmp = xs1
    xs1 = xs2
    xs2 = tmp
  }
  n1 = length(xs1)
  n2 = length(xs2)
  
  ys = c(xs1, xs2)
  
  # 如果一个样本中有超过五个值位于另一个样本的范围之外，返回1。这意味着样本之间的方差不同。
  count5 = function (ys1, ys2) {
    stopifnot(length(ys1) == length(ys2))
    extr1 = sum((ys1 < min(ys2))) + sum((ys1 > max(ys2)))
    extr2 = sum((ys2 < min(ys1))) + sum((ys2 > max(ys1)))
    
    out = max(extr1, extr2)
    return(as.integer(out > 5))
  }
  
  # 运行排列检验。
  # 如果H0成立，对于ys中的“大多数”（取决于所需的p值）等大小的2个样本集，count 5必须返回0。
  rep = 1000
  n.min = n1
  n = length(ys)
  c5s = numeric(rep)
  size = floor(n / 2)
  
  for (i in 1:rep) {
    k = sample(1:n, size, replace = FALSE)
    ys1 = ys[k]
    ys2 = ys[-k]
    # 计算count 5测试返回1（即不同方差）的次数
    c5s[i] = count5(ys1, ys2)
  }
  
  hist(c5s, probability = TRUE)
  
  # 根据count 5测试的不同方差频率返回结果
  return(mean(c5s))
}

mean = 0
# 使用不同标准差的不同大小的随机正态分布样本进行算法测试
result1 <- test_equal_variance(rnorm(n1, mean, 1), rnorm(n2, mean, 1))
result2 <- test_equal_variance(rnorm(n1, mean, 10), rnorm(n2, mean, 1))
result3 <- test_equal_variance(rnorm(n1, mean, 100), rnorm(n2, mean, 1))
result1
result2
result3

## ----message=FALSE------------------------------------------------------------
  # 导入所需的包
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
# 第一步：设计函数，计算a的值
calculate_a <- function(N, b1, b2, b3, f0) {
  # 生成随机变量 X1, X2, X3
  set.seed(123)  # 设置随机种子，以确保可重现性
  X1 <- rpois(N, 1)
  X2 <- rexp(N, 1)
  X3 <- rbinom(N, 1, 0.5)
  
  # 计算 Y=1 的概率
  a <- uniroot(function(a) {
    P_Y1 <- exp(a + b1*X1 + b2*X2 + b3*X3) / (1 + exp(a + b1*X1 + b2*X2 + b3*X3))
    mean(P_Y1) - f0
  }, interval = c(-100, 100))$root
  
  return(a)
}

# 第二步：输入参数
N <- 1e6
b1 <- 0
b2 <- 1
b3 <- -1
f0 <- c(0.1, 0.01, 0.001, 0.0001)

# 计算 a 的值
a_values <- sapply(f0, function(f) calculate_a(N, b1, b2, b3, f))

# 第三步：绘制图形
df <- data.frame(-log(f0), a_values)
colnames(df) <- c("-log f0", "a")
ggplot(df, aes(x = a, y = -log(f0))) +
  geom_line() +
  geom_point() +
  labs(x = "a", y = "-log f0") +
  theme_minimal()

## -----------------------------------------------------------------------------
sds <- c(0.5, 1, 2, 4)
is <- 5000:5500

mc.laplace <- function(sd) {
  N <- 10000

  # 定义标准拉普拉斯分布函数
  df <- function(x) {
    1/2 * exp(-abs(x))
  }

  # 定义生成随机数的函数
  rg <- function(mean) {
    rnorm(n = 1, mean = mean, sd = sd)
  }

  # Metropolis-Hastings算法
  rw <- function(N, df, rg) {
    x <- numeric(N)
    x[1] <- rg(1)
    k <- 0
    us <- runif(N)

    for (i in 2:N) {
      xt <- x[i-1]
      y <- rg(xt)
      res <- df(y) / df(xt)
      if (us[i] <= res) {
        x[i] <- y
      } else {
        x[i] <- xt
        k <- k + 1
      }
    }
    print(k)

    return(x)
  }

  return(rw(N, df, rg))
}

# 绘制图形
xs <- mc.laplace(sd = sds[1])
par(mfrow = c(length(sds), 1))
par(mar = c(2, 4, 2, 1))  # 调整边距值以减小边距的大小
plot(is, xs[is], type = "l", ylab = paste("标准差 = ", sds[1], sep = ''))
for (i in 2:length(sds)) {
  xs <- mc.laplace(sd = sds[i])
  plot(is, xs[is], type = "l", ylab = paste("标准差 = ", sds[i], sep = ''))
}
par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
m <- 5000
burn <- 1000

x <- matrix(0, m, 2)

rho <- 0.9
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2) * sigma1
s2 <- sqrt(1-rho^2) * sigma2

mean12 <- function(x2) mu1 + rho * sigma1 / sigma2 * (x2 - mu2)
mean21 <- function(x1) mu2 + rho * sigma2 / sigma1 * (x1 - mu1)

x[1,] <- c(mu1, mu2)

for (i in 2:m) {
  xt <- x[i-1,]
  xt[1] <- rnorm(1, mean12(xt[2]), s1)
  xt[2] <- rnorm(1, mean21(xt[1]), s2)
  x[i,] <- xt
}

x <- x[(burn+1):m,]

x <- data.frame(x)
lin.reg <- lm(X1 ~ X2, data = x)

par(mfrow = c(1, 2))
plot(x, cex = 0.5, main = "生成的数据")  
hist(lin.reg$residuals, main = "线性模型的残差")  
par(mfrow = c(1, 1))

## ----message=FALSE------------------------------------------------------------
  library(coda)

## -----------------------------------------------------------------------------
l <- 10
s <- c(1, 2, 3, 4)
k <- length(s)

# 定义Rayleigh分布函数
sigma <- 4
df <- function(y, sigma) {
  if (any(y < 0)) return(0)
  stopifnot(sigma > 0)
  return(y / sigma^2 * exp(-y^2 / (2 * sigma^2)))
}

dg <- function(y, xt) {
  dchisq(y, df = xt)
}

rg <- function(xt) {
  rchisq(1, df = xt)
}

# Metropolis-Hastings算法
mh <- function(s, l) {
  x <- numeric(l)
  us <- runif(l)
  x[1] <- s
  for (i in 2:l) {
    xt <- x[i - 1]
    y <- rg(xt)
    res <- df(y, sigma) / df(xt, sigma) * dg(xt, y) / dg(y, xt)
    if (us[i] <= res) {
      x[i] <- y
    } else {
      x[i] <- xt
    }
  }
  return(x)
}

xs <- matrix(sapply(1:k, function(i) mh(s[i], l)), nrow = k, byrow = TRUE)

# 可视化生成的样本
plotHists <- function() {
  par(mfrow = c(1, k))
  for (i in 1:k) {
    hist(xs[i, ], probability = TRUE, breaks = 100)
    x.axis <- seq(min(xs[i, ]), max(xs[i, ]), by = 0.01)
    lines(x.axis, df(x.axis, sigma))
  }
  par(mfrow = c(1, 1))
}

plotChains <- function() {
  burn <- 2000
  is <- (burn + 1):l
  par(mfrow = c(1, k))
  for (i in 1:k) {
    plot(is, xs[i, is], type = "l")
  }
  par(mfrow = c(1, 1))
}

# Gelman-Rubin诊断方法
gelman.rubin <- function(psis) {
  psi.means <- rowMeans(psis)
  B <- l * var(psi.means)
  W <- mean(apply(psis, MARGIN = 1, "var"))
  var.hat <- (l - 1) / l * W + 1 / l * B
  return(var.hat / W)
}

div <- matrix(sapply(1:k, function(i) 1:l), nrow = k, byrow = TRUE)
psis <- t(apply(xs, MARGIN = 1, "cumsum")) / div

r.hats <- sapply(2:l, function(j) gelman.rubin(psis[, 1:j]))

# 绘制Gelman-Rubin诊断结果
plot(2:l, r.hats, type = "l")
abline(h = 1.2)

## -----------------------------------------------------------------------------
# 定义似然函数
likelihood <- function(lambda, data) {
  likelihood <- prod(exp(-lambda * data[, 1]) - exp(-lambda * data[, 2]))
  return(-likelihood)  # 求最大值，因此取负号
}

# 使用直接极大化观测数据的似然函数求解λ的MLE
direct_mle <- function(data) {
  result <- optim(par = 1.0, fn = likelihood, data = data, method = "BFGS")  # 最大化似然函数
  lambda_mle <- result$par  # 得到λ的MLE
  return(lambda_mle)
}

# 使用EM算法求解λ的MLE
em_algorithm <- function(data) {
  lambda_initial <- 1.0  # 初始参数值
  result <- optim(par = lambda_initial, fn = likelihood, data = data, method = "BFGS")  # 最大化似然函数
  lambda_mle <- result$par  # 得到λ的MLE
  return(lambda_mle)
}

# 给定观测数据
data <- matrix(c(11, 12, 8, 9, 27, 28, 3, 14, 6, 17, 0, 1, 23, 24, 10, 11, 24, 25, 2, 3), ncol = 2, byrow = TRUE)

# 直接极大化观测数据的似然函数
lambda_mle_direct <- direct_mle(data)
cat("直接极大化观测数据的似然函数得到的λ的MLE：", lambda_mle_direct, "\n")

# 采用EM算法求解λ的MLE
lambda_mle_em <- em_algorithm(data)
cat("采用EM算法得到的λ的MLE：", lambda_mle_em, "\n")

## -----------------------------------------------------------------------------
dim(c(1, 2, 3))

dim(c('a', 'b', 'c'))

## -----------------------------------------------------------------------------
#请注意，矩阵是原子的（不能包含不同的类型）。因此我们会期望进行类型强制转换。
as.matrix(data.frame(a = 1, b = 'char'))

#注意，1 已经变成了一个字符型，"1"。

## -----------------------------------------------------------------------------
df <- data.frame()
df

class(df)

dim(df)

#能否有0行但是>=1列的数据框
df <- data.frame(1, seq(10))[FALSE, ]
dim(df)

#或者有>=1行但是0列的数据框
df <- data.frame(1, seq(10))[, FALSE]
dim(df)


## -----------------------------------------------------------------------------
### The function below scales a vector so it falls in the range [0, 1]. 

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

### How would you apply it to every column of a data frame?

lapply(mtcars, scale01)
# lapply goes by column on a data frame by default.

### How would you apply it to every numeric column in a data frame?

# for a data set like iris, we get an error when applying because of the non-numeric data frames.
iris[sapply(iris, is.numeric)]  # select only numeric cols.
lapply(iris[sapply(iris, is.numeric)], scale01)  # apply a function to only the numeric cols.

## -----------------------------------------------------------------------------
# a

numeric_iris <- iris[-5]
vapply(numeric_iris, sd, numeric(1))

# b

df_sd <- function(df) vapply(df[vapply(df, is.numeric, logical(1))], sd, numeric(1))
df_sd(iris)
df_sd(mtcars)

## ----warning=FALSE------------------------------------------------------------
n <- 100
a <- 30
b <- 60

df <- function(x, y) {
  # 一般二项式系数
  choose(n, x) * y^(x + a - 1) * (1 - y)^(n - x + b - 1)
}

m <- 10000
d <- 2

x <- matrix(0, nrow = m, ncol = d)

for (i in 2:m) {
  xt <- x[i - 1, ]
  xt[1] <- rbinom(1, n, xt[2]) # 从二项分布生成 x 的值
  xt[2] <- rbeta(1, xt[1] + a, n - xt[1] + b) # 从贝塔分布生成 y 的值
  x[i, ] <- xt
}

plot(x, cex = 0.1)
xs <- seq(min(x[, 1]), max(x[, 1]), length.out = 200)
ys <- seq(min(x[, 2]), max(x[, 2]), length.out = 200)
zs <- outer(xs, ys, Vectorize(df))

# 绘制密度轮廓图
contour(xs, ys, zs, add = TRUE, col = 2)

## -----------------------------------------------------------------------------
# R函数
myFunctionR <- function(n, a, b, m) {
  df <- function(x, y) {
    choose(n, x) * y^(x + a - 1) * (1 - y)^(n - x + b - 1)
  }
  
  x <- matrix(0, nrow = m, ncol = 2)
  
  for (i in 2:m) {
    xt <- x[i - 1, ]
    xt[1] <- rbinom(1, n, xt[2])
    xt[2] <- rbeta(1, xt[1] + a, n - xt[1] + b)
    x[i, ] <- xt
  }
  
  return(x)
}

## -----------------------------------------------------------------------------
#// Rcpp函数
#include <Rcpp.h>
#using namespace Rcpp;

#// [[Rcpp::export]]
#NumericMatrix myFunctionCpp(int n, int a, int b, int m) {
#  Function choose("choose");
#  NumericMatrix x(m, 2);
  
#  for (int i = 1; i < m; i++) {
#    NumericVector xt = x(i - 1, _);
#    xt[0] = as<int>(rbinom(1, n, xt[1]));
#    xt[1] = as<double>(rbeta(1, xt[0] + a, n - xt[0] + b));
#    x(i, _) = xt;
#  }
  
#  return x;
#}

## -----------------------------------------------------------------------------
library(microbenchmark)

n <- 100
a <- 30
b <- 60
m <- 10000

# 测试R函数的计算时间
benchmarkR <- microbenchmark(myFunctionR(n, a, b, m))

# 测试Rcpp函数的计算时间
#benchmarkCpp <- microbenchmark(myFunctionCpp(n, a, b, m))

# 打印计算时间结果
print(benchmarkR)
#print(benchmarkCpp)

