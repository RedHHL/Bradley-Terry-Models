
library(dplyr)

soccer <- read.csv('premier-league-matches.csv')
soccer$Date <- as.Date(soccer$Date, "%Y-%m-%d")
team_appearances <- soccer %>% 
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2023) %>% 
  group_by(Home) %>% 
  summarise(Distinct_Seasons = n_distinct(Season_End_Year)) %>% 
  filter(Distinct_Seasons >= 10)

soccer_filtered <- soccer %>% 
  filter(Home %in% team_appearances$Home & Away %in% team_appearances$Home)

######################################################################
calculate_rbar <- function(data, year, home_or_away) {
  previous_season_matches <- data %>% 
  filter(Season_End_Year == year)
  if (home_or_away == "home") {
    rbar <- mean(ifelse(previous_season_matches$FTR == "H", 3, ifelse(previous_season_matches$FTR == "D", 1, 0)))
  } else {
    rbar <- mean(ifelse(previous_season_matches$FTR == "A", 3, ifelse(previous_season_matches$FTR == "D", 1, 0)))
  }
}


calculate_team_ability <- function(beta, lambda, r_bar, past_results) {
  if (length(past_results) == 0) {
    return(beta * r_bar)
  }
  ability <- beta * r_bar * ((1 - lambda)^(length(past_results)))
  for (k in 1:length(past_results)) {
    ability <- ability + beta * lambda * ((1 - lambda)^(k-1)) * past_results[length(past_results) - k + 1]
  }
  return(ability)
}


get_past_results <- function(data, team, date, home_or_away) {
  if (home_or_away == "home") {
    past_matches <- data[data$Home == team & data$Date < date, ]
    past_results <- ifelse(past_matches$FTR == "H", 3, ifelse(past_matches$FTR == "D", 1, 0))
  } else {
    past_matches <- data[data$Away == team & data$Date < date, ]
    past_results <- ifelse(past_matches$FTR == "A", 3, ifelse(past_matches$FTR == "D", 1, 0))
  }

  return(past_results)
}

log_likelihood <- function(params, data, data2, year) {
  lambda1 <- params[1]
  lambda2 <- params[2]
  beta1 <- params[3]
  beta2 <- params[4]
  delta <- params[5]

  log_like <- 0
  rh_bar <- calculate_rbar(data2, year, "home")
  rv_bar <- calculate_rbar(data2, year, "away")

  
  for (i in 1:nrow(data)) {
    match <- data[i, ]
    past_home_results <- get_past_results(data, match$Home, match$Date, "home")
    past_away_results <- get_past_results(data, match$Away, match$Date, "away")
    
    home_ability <- calculate_team_ability(beta1, lambda1, rh_bar, past_home_results)
    away_ability <- calculate_team_ability(beta2, lambda2, rv_bar, past_away_results)

    
    logit_D <- delta + home_ability - away_ability
    logit_A <- -delta + home_ability - away_ability
    #print(logit_D)
    
    cum_prob_D <- exp(logit_D) / (1 + exp(logit_D))
    cum_prob_A <- exp(logit_A) / (1 + exp(logit_A))
    #print(cum_prob_A)
    
    prob_H <- 1 - cum_prob_D
    prob_D <- cum_prob_D - cum_prob_A
    prob_A <- cum_prob_A

    
    if (match$FTR == "H") {
      log_like <- log_like + log(prob_H)
    } else if (match$FTR == "D") {
      log_like <- log_like + log(prob_D)
    } else if (match$FTR == "A") {
      log_like <- log_like + log(prob_A)
    }
  }
  
  return(-log_like)
}



pre_covid_train  <- soccer_filtered %>%
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2018)
covid_train  <- soccer_filtered %>%
  filter(Season_End_Year >= 2020 & Season_End_Year <= 2021)
pre_covid_test  <- soccer_filtered %>%
  filter(Season_End_Year == 2019)
covid_test  <- soccer_filtered %>%
  filter(Season_End_Year == 2022)

pre_covid_train_fit <- nlminb(c(lambda1=0.003, lambda2=0.0014, beta1=-16, beta2=-26, delta=0.56),objective = log_likelihood, data = pre_covid_train, data2 = soccer_filtered, year = 2013, lower = c(0.00000000000001, 0.00000000000001, -Inf, -Inf,-Inf) ,upper = c(0.999999999999, 0.999999999999, Inf,Inf,Inf), control=list(trace=1,x.tol=0.0000000001))
covid_train_fit <- nlminb(c(lambda1=0.05, lambda2=0.05, beta1=0.65, beta2=1.18, delta=0.512),objective = log_likelihood, data = covid_train, data2 = soccer_filtered, year = 2019, lower = c(0.00000000000001, 0.00000000000001, -Inf, -Inf,-Inf) ,upper = c(0.999999999999, 0.999999999999, Inf,Inf,Inf), control=list(trace=1,x.tol=0.0000000001))


print(pre_covid_train_fit$par)
print(covid_train_fit$par)
######################################################################
library(ggplot2)

calculate_abilities_for_team <- function(team, soccer_data, data, mle_params, year) {
  ability_data <- data.frame(Date = as.Date(character()), HomeAbility = numeric(), AwayAbility = numeric())
  rh_bar <- calculate_rbar(data, year, "home")
  rv_bar <- calculate_rbar(data, year, "away")
  
  
  for (date in unique(soccer_data$Date)) {
    past_home_results <- get_past_results(soccer_data, team, date, "home")
    past_away_results <- get_past_results(soccer_data, team, date, "away")
    
    home_ability <- -calculate_team_ability(mle_params["beta1"], mle_params["lambda1"], rh_bar, past_home_results)
    away_ability <- -calculate_team_ability(mle_params["beta2"], mle_params["lambda2"], rv_bar, past_away_results)
    
    ability_data <- rbind(ability_data, data.frame(Date = date, HomeAbility = home_ability, AwayAbility = away_ability))
  }
  
  ability_data <- ability_data[order(ability_data$Date), ]
  ability_data$Date <- as.Date(ability_data$Date)
  return(ability_data)
}


plot_abilities <- function(abilities_data, team, plot_label) {
  p <-ggplot(abilities_data, aes(x = Date)) +
    geom_line(aes(y = HomeAbility), linetype="solid") + 
    geom_line(aes(y = AwayAbility), linetype="dotted") +  
    labs(title=plot_label,y = "Ability", x = "Date") +
    theme_minimal() +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none",
          legend.key = element_blank())
  
  return(p)
}
labels <- paste0(LETTERS[1:10])


library(patchwork)


pre_covid_team_plots <- list()
covid_team_plots <- list()
teams <- unique(soccer_filtered$Home)
label_index <- 1
for (team in teams) {
  pre_covid_team_data <- calculate_abilities_for_team(team, pre_covid_train, soccer_filtered, pre_covid_train_fit$par, 2013)
  pre_covid_team_plot <- plot_abilities(pre_covid_team_data, team, labels[label_index])
  pre_covid_team_plots[[label_index]]  <- pre_covid_team_plot
  covid_team_data <- calculate_abilities_for_team(team, covid_train, soccer_filtered, covid_train_fit$par, 2019)
  covid_team_plot <- plot_abilities(covid_team_data, team, labels[label_index])
  covid_team_plots[[label_index]]  <- covid_team_plot
  label_index <- label_index + 1
}

pre_covid_combined_plot <- wrap_plots(pre_covid_team_plots, ncol = 3)

covid_combined_plot <- wrap_plots(covid_team_plots, ncol = 3)

######################################################################


calculate_probabilities <- function(data, params, data1, year, get_past_results, calculate_team_ability) {
  predictions <- data.frame(A = numeric(nrow(data)), D = numeric(nrow(data)), H = numeric(nrow(data)))
  rh_bar <- calculate_rbar(data1, year, "home")
  rv_bar <- calculate_rbar(data1, year, "away")
  
  for (i in 1:nrow(data)) {
    match <- data[i, ]
    home_ability <- calculate_team_ability(params["beta1"], params["lambda1"], rh_bar, get_past_results(data, match$Home, match$Date, "home"))
    away_ability <- calculate_team_ability(params["beta2"], params["lambda2"], rv_bar, get_past_results(data, match$Away, match$Date, "away"))
    
    logit_H <- params["delta"] + home_ability - away_ability
    logit_A <- -params["delta"] + home_ability - away_ability
    
    prob_H <- 1 - exp(logit_H) / (1 + exp(logit_H))
    prob_A <- exp(logit_A) / (1 + exp(logit_A))
    prob_D <- 1 - prob_H - prob_A
    
    predictions[i, ] <- c(A = prob_A, D = prob_D, H = prob_H)
  }
  return(predictions)
}

calculate_scores <- function(data, params, data1, year, get_past_results, calculate_team_ability, log_likelihood) {
  predicted_probs <- calculate_probabilities(data, params, data1, year, get_past_results, calculate_team_ability)
  # Brier Score
  outcomes <- c("A", "D", "H")
  actual_probs <- matrix(as.numeric(unlist(lapply(data$FTR, function(x) outcomes == x))), nrow = nrow(data), byrow = TRUE)
  brier_score <- mean(rowSums((predicted_probs - actual_probs) ^ 2) / ncol(actual_probs))

  # Ranked Probability Score (RPS)
  predicted_cum_probs <- t(apply(predicted_probs, 1, cumsum))
  actual_cum_probs <- t(apply(actual_probs, 1, cumsum))
  rps <- mean(rowSums((predicted_cum_probs - actual_cum_probs) ^ 2) / (ncol(actual_probs) - 1))

  # AIC
  ll <- log_likelihood(params, data, data1, year)
  k <- length(params)
  aic <- 2 * k + 2 * ll
  bic <- log(length(data)) * k + 2 * ll

  return(list(BrierScore = brier_score, RankedProbabilityScore = rps, AIC = aic, BIC = bic,loglikehood = ll))
}

pre_covid_scores <- calculate_scores(
  data = pre_covid_test,
  params = pre_covid_train_fit$par,
  data1 = soccer_filtered,
  year = 2018,
  get_past_results = get_past_results,
  calculate_team_ability = calculate_team_ability,
  log_likelihood = log_likelihood
)

covid_scores <- calculate_scores(
  data = covid_test,
  params = covid_train_fit$par,
  data1 = soccer_filtered,
  year = 2021,
  get_past_results = get_past_results,
  calculate_team_ability = calculate_team_ability,
  log_likelihood = log_likelihood
)


print(pre_covid_scores)
print(covid_scores)