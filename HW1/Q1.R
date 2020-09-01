library(ggplot2)
library(RColorBrewer)
library(scales)

N = 20000

rf <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))
r <- rf(32)

plot_joint = function(x, y) {
  df = data.frame(x, y)
  p = ggplot(df, aes(x=x, y=y)) + stat_density_2d(aes(fill=..level..), geom="polygon") + scale_fill_gradientn(colours=r) + ggtitle("joint")
  print(p)
}

MARGIN = 1e-1
plot_conditional = function(first, second, horizontal_axis_title, title_template) {
  df = data.frame(first, second)
  for (val in c(0, 1, 2)) {
    title = paste(title_template, val, sep=" ")
    new_df = subset(df, second > val - MARGIN & second < val + MARGIN, select=first)
    p = ggplot(new_df, aes(x=first)) + geom_histogram(aes(y=..density..), binwidth=0.01) + labs(x=horizontal_axis_title) + ggtitle(title)
    print(p)
  }
}

plot_intervention_on_y = function(x) {
  sampled_x = data.frame(sample(x, size=10000))
  names(sampled_x)[1] = 'x'
  p = ggplot(sampled_x, aes(x=x)) + geom_histogram(aes(y=..density..), binwidth=0.01) + labs(x='x') +  ggtitle("interved on y = anything")
  print(p)
}

plot_everything = function(x, y) {
  plot_joint(x, y)  # p(x, y)
  plot_conditional(y, x, horizontal_axis_title='y', title_template='given x =')  # p(y|x)
  plot_conditional(x, y, horizontal_axis_title='x', title_template='given y =')  # p(x|y)
  plot_conditional(y, x, horizontal_axis_title='y', title_template='intervened on x =')  # p(y, do(x)), intervention on cause is like conditioning on cause
  plot_intervention_on_y(x)  # p(x, do(y))
}

# Main
x = rnorm(N, mean=0, sd=1)

# Part 1
n_y = rnorm(N, mean=0, sd=1)
y = x^3 + n_y
plot_everything(x, y)

# Part 2
n_y = rt(N, df=1)
y = 2*x + n_y
plot_everything(x, y)

# Part 3
n_y = rt(N, df=5)
y = 2*x + n_y
plot_everything(x, y)

# Part 4
n_y = rt(N, df=20)
y = 2*x + n_y
plot_everything(x, y)
