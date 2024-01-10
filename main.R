source('algorithm.r')
options(scipen = 999)

simulate <- function(delta, transition_matrix, lambda, size){
  states <- 1: nrow(transition_matrix)
  current_state <- sample(states, 1, prob = delta)
  simulated_data <- integer(size)
  
  for (i in 1:size) {
    simulated_data[i] <- rpois(1, lambda[current_state])
    next_possible_states_probs <- transition_matrix[current_state, ]
    current_state <- sample(states, 1, prob = next_possible_states_probs)
  }
  
  return(simulated_data)
  
}

random_probability_vector <- function(dimension){
  unormalized_vector<- runif(dimension, 0, 1)
  return(unormalized_vector / sum(unormalized_vector))
  
}

random_transition_matrix <- function(dimension){
  transition_matrix <- matrix(0, nrow = dimension, ncol = dimension)
  
  for (i in 1:dimension) {
    transition_matrix[i,] <- random_probability_vector(dimension)
  }
  return(transition_matrix)
}


data <- read.table('earthquakes.txt')
data

colnames(data) <- c('date', 'earth_quakes_count')


observations <- data$earth_quakes_count
mean(observations)
var(observations)

hist(observations, probability = TRUE, ylim = c(0, .10))
lines(0:max(observations), dpois(0:max(observations), mean(observations)), col = 'red')

acf(observations)

observations[2:107]


transition_matrix <- random_transition_matrix(6)
transition_matrix
delta <- random_probability_vector(6)
lambda <- c(20, 10, 15, 17, 18, 5)


estimations <- em.poisson.hmm(delta, transition_matrix, observations, lambda, n_iter = 500)


simulated_data <- simulate(as.vector(estimations$delta), estimations$transition_matrix, estimations$lambda, 107)

hist(simulated_data, probability = TRUE)

hist(observations)


simulated_data_ecdf <- ecdf(simulated_data)
observations_ecdf <- ecdf(observations)
plot(observations_ecdf, verticals=TRUE, do.points=FALSE, col='blue')
plot(simulated_data_ecdf, verticals=TRUE, do.points=FALSE, add=TRUE, col='red')

summary(observations)
summary(simulated_data)

