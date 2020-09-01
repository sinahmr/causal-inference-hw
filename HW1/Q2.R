#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("graph", "RBGL"))

library(kpcalg)
library(mgcv)

N = 200

print_dir = function(x, y) {
  df = data.frame(x, y)
  
  model = gam(y ~ s(x), data=df)  # Non-linear regression (automatically finds the “best” knots)
  predictions = predict(model, df)
  hsic.perm(x, y - predictions)
  p_val_1 = hsic.perm(x, y - predictions)$p.value
  
  model = gam(x ~ s(y), data=df)
  predictions = predict(model, df)
  hsic.perm(y, x - predictions)
  p_val_2 = hsic.perm(y, x - predictions)$p.value
  
  if (p_val_1 > p_val_2) {
    print("X -> Y")
  } else {
    print("Y -> X")
  }
}

model_a = function() {
  print('Model A:')
  for (i in 1:10) {
    x = rnorm(N, mean=0, sd=1)
    n_y = rnorm(N, mean=0, sd=1)
    y = x^3 + n_y
    
    print_dir(x, y)
  }
  print('---')
}

model_b = function(v) {  # v is the parameter of the t-dist
  print(paste('Model B, v =', v, ':', sep=" "))
  for (i in 1:10) {
    x = rnorm(N, mean=0, sd=1)
    n_y = rt(N, df=v)
    y = 2*x + n_y
    
    print_dir(x, y)
  }
  print('---')
}

tuebingen = function() {
  df = read.delim("Causal Inference/HW1/pair0001.txt", header=FALSE, sep=" ", dec = ".")
  print_dir(df$V1, df$V2)
  
  df = read.delim("Causal Inference/HW1/pair0002.txt", header=FALSE, sep=" ", dec = ".")
  print_dir(df$V1, df$V2)
}
model_a()
model_b(1)
model_b(5)
model_b(20)
tuebingen()
