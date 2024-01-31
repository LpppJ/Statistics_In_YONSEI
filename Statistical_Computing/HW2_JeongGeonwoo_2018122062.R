# STA3127 Statistical Computing
# HW2
# Name : Jeong Geonwoo
# ID : 2018122062

rm(list=ls())
set.seed(2018122062)

##### Q1 #####
u_grid <- seq(0, 1, length.out=1001)
v_grid <- seq(0, 1, length.out=1001)

integrand1 <- function(u,v){
  value <- 2*sin(20*u^2+8*u*v+v^2) * 1*(0<=u & u<=1) * 1*(0<=v & v<=1)
  return (value)
}

z_values1 <- matrix(NA, nrow = length(u_grid), ncol = length(v_grid))
for (i in 1:length(u_grid)) {
  for (j in 1:length(v_grid)) {
    z_values1[i, j] <- integrand1(u_grid[i], v_grid[j])
  }
}

contour(u_grid, v_grid, z_values1, xlab = "u", ylab = "v")



integrand2 <- function(u,v){
  value <- 10*sin(4*u^2+25*v^2) * 1*(0<=u & u<=1) * 1*(4/5*u<=v & v<=4/5*u+1/5)
  return (value)
}
z_values2 <- matrix(NA, nrow = length(u_grid), ncol = length(v_grid))
for (i in 1:length(u_grid)) {
  for (j in 1:length(v_grid)) {
    z_values2[i, j] <- integrand2(u_grid[i], v_grid[j])
  }
}
contour(u_grid, v_grid, z_values2, xlab = "u", ylab = "v")

integral1_list <- c()
integral2_list <- c()
for (i in c(1:1000)){

  integral1 = 0
  integral2 = 0
  for (i in 1:100) {
    u <- runif(1, 0, 1)
    v <- runif(1, 0, 1)
    integral1 = integral1 + integrand1(u, v)
    integral2 = integral2 + integrand2(u, v)
  }
  integral1 = integral1 / 100
  integral2 = integral2 / 100
  integral1_list <- c(integral1_list, integral1)
  integral2_list <- c(integral2_list, integral2)
}

cat("integral1 mean:", mean(integral1_list), ", integral1 variance:", var(integral1_list))
cat("integral2 mean:", mean(integral2_list), ", integral2 variance:", var(integral2_list))



##### Q2 #####
prob_no_bomb<- function(p, q, k, n){
  ########################################
  # input
  # p : the number of rows of map
  # q : the number of columns of map
  # k : the number of bombs in the map
  # n : the number of realizations of the indicator variable
  
  # output
  # success_ratio : estimated expectation(probability) of "zero bombs in the vicinity of the selected cell"
  ########################################
  
  count_success = 0
  # if clicked cell and surrounding 8 cells have no bomb,
  # then "count_success" increases by 1
  
  for (realization in c(1:n)){
    # Empty map
    map <- matrix(0, nrow = p, ncol = q)
    
    # Location of bombs : Initializing
    bomb_location <- c(1:(p*q))
    # Location of bombs : Random Permutation
    for (index in c(1:k)){
      I = floor(p*q*runif(1,0,1))+1
      temp <- bomb_location[I]
      bomb_location[I] <- bomb_location[index]
      bomb_location[index] <- temp
    }
    # Location of bombs : Select k-location
    bomb_location = bomb_location[1:k]
    
    # Placing the bombs
    for (location in bomb_location) {
      row_idx <- (location - 1) %/% q + 1
      col_idx <- (location - 1) %% q + 1
      map[row_idx, col_idx] <- 1
    }
    # Initial click
    click <- floor(p*q*runif(1,0,1))+1
    click_p <- (click-1) %/% q + 1
    click_q <- (click-1) %% q + 1
    
    # Count the number of bombs in surrounding cells
    count_bombs = 0
    for (i in c(-1,0,1)){
      for (j in c(-1,0,1)){
        if ((1 <= click_p+i & click_p+i <= p) & (1<= click_q+j & click_q+j <= q)){
          count_bombs = count_bombs + map[click_p+i, click_q+j]
        }
      }
    }
    # if there is no bomb in 9-cell, count_success increases by 1
    if (count_bombs == 0){
      count_success = count_success + 1
    }
  }
  success_ratio = count_success / n
  return (success_ratio)
}

prob_no_bomb(1,1,1,1000)
prob_no_bomb(3,3,1,1000)
prob_no_bomb(20,15,40,1000)
