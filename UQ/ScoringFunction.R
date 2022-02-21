# Uncertainty Quantification: Scoring function
# Federica Milinanni (fedmil@kth.se), based on:
# Copyright (C) 2021 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Joao Antunes (joaodgantunes@gmail.com) and Olivia Eriksson (olivia@kth.se)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

getScoreTimeSeries <- function(yy_sim, yy_exp){
  #if one dimensional output
  if(class(yy_sim) == "numeric"){
    #if one dimensional output
    distance <- mean((yy_sim-yy_exp)^2)
  } else {
    #if multidimensional output
    p <- 2
    q <- 2
    weights_vec <- rep(1,dim(yy_sim)[1])
    distance <- sum(((colSums(((yy_sim-yy_exp)*weights_vec)^p))^(1/p))^q)/dim(yy_sim)[2]
  }
  return(distance)
}


getScoreTimeSeriesErrorNormalisation <- function(yy_sim, yy_exp, err){
  #if one dimensional output
  if(class(yy_sim) == "numeric"){
    #if one dimensional output
    distance <- sum(((yy_sim-yy_exp)/err)^2)/length(yy_sim)
  } else {
    #if multidimensional output
    p <- 2
    q <- 2
    weights_vec <- rep(1,dim(yy_sim)[1])
    distance <- sum(((colSums((((yy_sim-yy_exp)/err)*weights_vec)^p))^(1/p))^q)/dim(yy_sim)[2]
  }
  return(distance)
}


getScore  <- function(ytarget, yy, ytarget_min, ytarget_max, ymin, ymax){
  distance <- 100
  yy_n <- (yy-ymin)/(ymax-ymin)
  ytarget_n <- (ytarget-ytarget_min)/(ytarget_max-ytarget_min)
  
  if (!(all(yy==0, na.rm=T) | any(yy<0, na.rm=T)| any(is.na(yy))| any(yy>120))){
    distance <- 0
    
    for (j in 1:length(ytarget)){
      sigma <- 1
      distance <- distance + ((yy_n[j]-ytarget_n[j])/sigma)^2
    }
    distance <- distance / length(ytarget_n)
  }

  distance
}
