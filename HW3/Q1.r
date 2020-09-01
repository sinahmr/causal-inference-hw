library("pcalg")
library("lmtest")

df = read.delim("./model02.txt", header=FALSE, sep=",", dec = ".")
colnames(df) = c("Y1", "Y2", "Y3")

# CASE1: TIME SERIES
grangertest(Y1 ~ Y2, order = 1, data = df)
grangertest(Y1 ~ Y3, order = 1, data = df)
grangertest(Y2 ~ Y1, order = 1, data = df)
grangertest(Y2 ~ Y3, order = 1, data = df)
grangertest(Y3 ~ Y1, order = 1, data = df)
grangertest(Y3 ~ Y2, order = 1, data = df)

# CASE2: NOT TIME SERIES
augmented_df = cbind(df[-nrow(df), ], df[-1, ])
augmented_V = c("Y1old", "Y2old", "Y3old", "Y1new", "Y2new", "Y3new")
colnames(augmented_df) = augmented_V

# Find DAG
n = nrow (augmented_df)
fixed_gaps = matrix(FALSE, nrow = 6, ncol = 6); fixed_gaps[1:3, 1:3] = TRUE; fixed_gaps[4:6, 4:6] = TRUE

dag = pc(suffStat = list(C = cor(augmented_df), n = n), indepTest = gaussCItest, alpha=0.01, labels = augmented_V, fixedGaps = fixed_gaps, verbose = FALSE)
summary(dag)

# dag = fci(list(C = cor(augmented_df), n = n), indepTest=gaussCItest, alpha = 0.01, labels = augmented_V, fixedGaps = fixed_gaps, doPdsep = FALSE)
# summary(dag)

# dag = lingam(augmented_df)
# as(dag, "amat")

# Find Coefficients
lm("Y2new ~ Y1old", data=augmented_df)
