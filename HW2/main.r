library(gtools)
library(pcalg)

# pval -> 1: less edge deletion
# alpha -> 1: less edge deletion

PART_A_CORRECT_ADJ = rbind(
  c(0, 0, 1, 0, 0, 0),
  c(0, 0, 1, 0, 0, 1),
  c(0, 0, 0, 1, 1, 0),
  c(0, 0, 0, 0, 1, 0),
  c(0, 0, 0, 0, 0, 1),
  c(0, 0, 0, 0, 0, 0)
)

my_combn = function(arr, l) {  # combn is stupid when arr is a single element
  if (length(arr) == 1) {
    return(list(arr))
  }
  return(combn(arr, l, simplify = TRUE))
}

is_independent = function(df, names, i, j, pval) {
  x = names[i]
  y = names[j]
  res = cor.test(df[[x]], df[[y]], method="pearson")
  if (res$p.value < pval)
    return(FALSE)
  return(TRUE)
}

partial_correlation = function(df, names, i, j, k) {  # k is a list
  x = names[i]
  y = names[j]
  z = names[k]
  rhs = paste(z, collapse=" + ")
  
  fmla_x = as.formula(paste(x, "~", rhs, sep=" "))
  res_x = lm(fmla_x, data=df)$residuals
  
  fmla_y = as.formula(paste(y, "~", rhs, sep=" "))
  res_y = lm(fmla_y, data=df)$residuals
  
  rho = cor.test(res_x, res_y, method="pearson")$estimate
  return(rho)
}

is_conditionally_independent = function(df, names, i, j, k, alpha) {
  n = nrow(df)
  rho = partial_correlation(df, names, i, j, k)
  z = 0.5 * (log(1+rho) - log(1-rho))
  threshold = qnorm(1 - alpha/2, mean=0, sd=1)
  
  if (sqrt(n -length(z) - 3) * abs(z) > threshold)
    return(FALSE)
  return(TRUE)
}

delete_edge <- defmacro(
  adj, S, i, j, k, expr={
    adj[i, j] = 0
    adj[j, i] = 0
    for (x in k) {
      S[i, j, x] = 1
      S[j, i, x] = 1
    }
    }
  )

PC = function(df, pval, alpha) {
  names = colnames(df)
  n_V = length(names)
  adj = array(1, dim=c(n_V, n_V))  # adj is matrix C of the algorithm
  diag(adj) = 0
  S = array(0, dim=c(n_V, n_V, n_V))
  
  for (l in 0:(n_V-2)) {
    for (i in 1:n_V) {
      neighbors = which(adj[i,] == 1)
      n_neighbors = length(neighbors)
      if (n_neighbors - 1 < l) 
        next
      
      for (idx in 1:n_neighbors) {
        j = neighbors[idx]
        if (l == 0) {
          if (is_independent(df, names, i, j, pval)) {
            delete_edge(adj, S, i, j, c())
          }
        } else {
          matrix_of_ks = my_combn(neighbors[-idx], l)  # a[-1] means every elements of a except for the first element
          for (t in ncol(matrix_of_ks)) {
            k = matrix_of_ks[, t]
            if (is_conditionally_independent(df, names, i, j, k, alpha)) {
              delete_edge(adj, S, i, j, k)
            }
          }
        }
      }
    }
  }
  return(list(adj=adj, S=S))
}


apply_v_structures = function(adj, S) {
  n_V = ncol(adj)
  for (x in 1:n_V) {
    for (z in (1:n_V)[-x]) {
      for (y in (1:n_V)[c(-x, -z)]) {
        if (adj[x, z]==1 && adj[z, x]==1 && adj[z, y]==1 && adj[y, z]==1 && adj[x, y]==0 && adj[y, x]==0) {
          if (S[x, y, z] == 0) {
            adj[z, x] = 0
            adj[z, y] = 0
          }
        }
      }
    }
  }
  return(adj)
}

meek = function(adj) {
  n_V = ncol(adj)
  flag = TRUE
  while (flag) {
    flag = FALSE
    for (x in 1:n_V) {
      for (y in (1:n_V)[-x]) {
        for (z in (1:n_V)[c(-x, -y)]) {
          if (adj[x, y]==1 && adj[y, x]==0 && adj[y, z]==1 && adj[z, y]==1 && adj[x, z]==0 && adj[z, x]==0) {
            # print(paste("type I", x, y, z, sep=" "))
            adj[z, y] = 0
            flag = TRUE
          }
          if (adj[x, y]==1 && adj[y, x]==0 && adj[y, z]==1 && adj[z, y]==1 && adj[z, x]==1 && adj[x, z]==0) {
            # print(paste("type II", x, y, z, sep=" "))
            adj[y, z] = 0
            flag = TRUE
          }
          for (w in (1:n_V)[c(-x, -y, -z)]) {
            if (adj[x, y]==1 && adj[y, x]==1 && adj[x, z]==1 && adj[z, x]==1 && adj[x, w]==1 && adj[w, x]==1 && adj[y, z]==0 && adj[z, y]==0 && adj[y, w]==1 && adj[w, y]==0 && adj[z, w]==1 && adj[w, z]==0) {
              # print(paste("type III", x, y, z, w, sep=" "))
              adj[w, x] = 0
              flag = TRUE
            }
            if (adj[x, w]==0 && adj[w, x]==0 && adj[x, y]==1 && adj[y, x]==1 && adj[x, z]==0 && adj[z, x]==1 && adj[y, z]==1 && adj[z, y]==1 && adj[y, w]==1 && adj[w, y]==1 && adj[z, w]==0 && adj[w, z]==1) {
              # print(paste("type IV", x, y, z, w, sep=" "))
              adj[x, y] = 0
              flag = TRUE
            }
          }
        }
      }
    }
  }
  return(adj)
}

cpdag_from_data = function(df, pval, alpha) {
  tmp = PC(df, pval, alpha)
  adj = tmp$adj
  S = tmp$S
  v_str_adj = apply_v_structures(adj, S)
  cpdag = meek(v_str_adj)
  return(cpdag)
}

dag_hamming = function(adj1, adj2) {
  return(sum(adj1 != adj2))
}

skeleton_hamming = function(adj1, adj2) {
  n_V = ncol(adj1)
  for (i in 1:n_V) {
    for (j in 1:n_V) {
      if (adj1[i, j] == 1) adj1[j, i] = 1
      if (adj2[i, j] == 1) adj2[j, i] = 1
    }
  }
  return(sum(adj1 != adj2))
}


# Part alef
N = 1000
eps1 = rnorm(N, mean=0, sd=1)
eps2 = rnorm(N, mean=0, sd=1)
eps3 = rnorm(N, mean=0, sd=1)
eps4 = rnorm(N, mean=0, sd=1)
eps5 = rnorm(N, mean=0, sd=1)
eps6 = rnorm(N, mean=0, sd=1)
x1 = 1.2*eps1
x2 = eps2
x3 = 2*x1 - 0.5*x2 + eps3
x4 = 0.4*x3 + eps4
x5 = 0.8*x3 - x4 + 0.6*eps5
x6 = x2 + x5 + eps6
df = data.frame(x1, x2, x3, x4, x5, x6)

cpdag = cpdag_from_data(df, 0.01, 0.01)


# Part be
n_V = 20
n = 1000
dag = randomDAG(n_V, 0.2, lB=0.1, uB=1, V=paste("x", seq(1:n_V), sep=""))
part_b_correct_adj = ifelse(as(dag, "matrix") > 0, 1, 0)  # ifelse is R's equivalent of numpy's where
df = as.data.frame(rmvDAG(n, dag, errDist="normal"))

alphas = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99)
dag_h = list()
skeleton_h = list()
i = 0
for (alpha in alphas) {
  i = i + 1
  cpdag = cpdag_from_data(df, 1e-7, alpha)
  dag_h[[i]] = dag_hamming(cpdag, part_b_correct_adj)
  skeleton_h[[i]] = skeleton_hamming(cpdag, part_b_correct_adj)
}

plot(alphas, dag_h, type="l", col="red", main="Hamming / Alpha", xlab="Alpha", ylab="Hamming")
lines(alphas, skeleton_h, col="blue")
legend("topleft", legend=c("DAG Hamming", "Skeleton Hamming"), col=c("red", "blue"), lty=1:1, cex=0.8)
