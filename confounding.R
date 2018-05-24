# P(y|do(x)) = P(y|x)
# x and y are not confounded whenever the observationally witnessed association between them is the same 
# as the association that would be measured in a controlled experiment, with x randomized

set.seed(1)
x = rnorm(100, 1, 3)
y = rnorm(100, 1, 3)

z = 2 * x + 3 * y + rnorm(100, 0, 3) 

summary(lm(y ~ x))
summary(lm(z ~ x))
summary(lm(z ~ y))
summary(lm(z ~ x + y))
