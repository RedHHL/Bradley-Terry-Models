library(dplyr)
library(flextable)
library(officer)

split_data_by_year <- function(data, train_start, train_end, test_start, test_end) {
  train_data <- data %>% filter(Season_End_Year >= train_start & Season_End_Year <= train_end)
  test_data <- data %>% filter(Season_End_Year >= test_start & Season_End_Year <= test_end)
  list(train = train_data, test = test_data)
}
soccer <- read.csv('premier-league-matches.csv')
soccer$FTR <- factor(soccer$FTR, levels = c("A", "D", "H"), ordered = TRUE)
soccer$Home <- factor(soccer$Home)
soccer$Away <- factor(soccer$Away)
reference_team <- "Liverpool"

team_appearances <- soccer %>% 
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2023) %>% 
  group_by(Home) %>% 
  summarise(Distinct_Seasons = n_distinct(Season_End_Year)) %>% 
  filter(Distinct_Seasons >= 10)
soccer_filtered <- soccer %>% 
  filter(Home %in% team_appearances$Home & Away %in% team_appearances$Home)

split_2014_2019 <- split_data_by_year(soccer_filtered, 2014, 2018, 2019, 2019)
split_2020_2022 <- split_data_by_year(soccer_filtered, 2020, 2021, 2022, 2022)
#############################################################################
log_likelihood <- function(params, data, reference, team) {
  abilities_home <- params[1:9]
  abilities_away <- params[10:18]
  delta <- params[19] 

  names(abilities_home) <- team
  names(abilities_away) <- team
  
  log_like <- 0
  
  for (i in 1:nrow(data)) {
    match <- data[i, ]
    
    home_ability <- ifelse(as.character(match$Home) == reference, 0, abilities_home[as.character(match$Home)])
    away_ability <- ifelse(as.character(match$Away) == reference, 0, abilities_away[as.character(match$Away)])
    
    logit_D <- delta + home_ability - away_ability
    logit_A <- -delta + home_ability - away_ability
    
    cum_prob_D <- exp(logit_D) / (1 + exp(logit_D))
    cum_prob_A <- exp(logit_A) / (1 + exp(logit_A))
    
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
teams <- setdiff(unique(c(soccer_filtered$Home, soccer_filtered$Away)), reference_team)
initial_params <- c(rep(1, 18), 1)
names(initial_params) <- c(paste0("Home", teams), paste0("Away", teams), "delta")
pre_covid_train_fit <- nlminb(initial_params,objective = log_likelihood, data = split_2014_2019$train, reference = reference_team, team = teams, control=list(trace=1,x.tol=0.0000000001))
covid_train_fit <- nlminb(initial_params,objective = log_likelihood, data = split_2020_2022$train, reference = reference_team, team = teams, control=list(trace=1,x.tol=0.0000000001))
post_covid_train_fit <- nlminb(initial_params,objective = log_likelihood, data = filter(soccer_filtered, Season_End_Year == 2023), reference = reference_team, team = teams, control=list(trace=1,x.tol=0.0000000001))
#############################################################################
calculate_probabilities <- function(data, params, team, reference) {
  predictions <- data.frame(A = numeric(nrow(data)), D = numeric(nrow(data)), H = numeric(nrow(data)))
  abilities_home <- params[1:9]
  abilities_away <- params[10:18]
  delta <- params[19] 
  names(abilities_home) <- team
  names(abilities_away) <- team


  for (i in 1:nrow(data)) {
    match <- data[i, ]
    home_ability <- ifelse(as.character(match$Home) == reference, 0, abilities_home[as.character(match$Home)])
    away_ability <- ifelse(as.character(match$Away) == reference, 0, abilities_away[as.character(match$Away)])

    
    
    logit_D <- delta + home_ability - away_ability
    logit_A <- -delta + home_ability - away_ability
    
    cum_prob_D <- exp(logit_D) / (1 + exp(logit_D))
    cum_prob_A <- exp(logit_A) / (1 + exp(logit_A))
    prob_H <- 1 - cum_prob_D
    prob_D <- cum_prob_D - cum_prob_A
    prob_A <- cum_prob_A
    
    
    predictions[i, ] <- c(A = prob_A, D = prob_D, H = prob_H)
  }
  return(predictions)
}

calculate_scores <- function(data, params, reference, team, log_likelihood) {
  predicted_probs <- calculate_probabilities(data, params, team, reference)
  # Brier Score
  outcomes <- c("A", "D", "H")
  actual_probs <- matrix(as.numeric(unlist(lapply(data$FTR, function(x) outcomes == x))), nrow = nrow(data), byrow = TRUE)
  brier_score <- mean(rowSums((predicted_probs - actual_probs) ^ 2) / ncol(actual_probs))
  
  # Ranked Probability Score (RPS)
  predicted_cum_probs <- t(apply(predicted_probs, 1, cumsum))
  actual_cum_probs <- t(apply(actual_probs, 1, cumsum))
  rps <- mean(rowSums((predicted_cum_probs - actual_cum_probs) ^ 2) / (ncol(actual_probs) - 1))
  
  # AIC
  ll <- log_likelihood(params, data, reference,team)
  
  k <- length(params)
  aic <- 2 * k + 2 * ll
  bic <- log(length(data)) * k + 2 * ll
  
  return(list(BrierScore = brier_score, RankedProbabilityScore = rps, AIC = aic, BIC = bic, loglikelihood=ll))
}

pre_covid_scores <- calculate_scores(
  data = split_2014_2019$test,
  params = pre_covid_train_fit$par,
  reference = reference_team,
  team = teams,
  log_likelihood = log_likelihood
)
covid_scores <- calculate_scores(
  data = split_2020_2022$test,
  params = covid_train_fit$par,
  reference = reference_team,
  team = teams,
  log_likelihood = log_likelihood
)
#############################################################################
extract_ranks <- function(estimate, reference_team) {
  home_abilities <- -estimate[grepl("Home", names(estimate))]
  away_abilities <- estimate[grepl("Away", names(estimate))]

  team_names <- gsub("Home|Away", "", names(home_abilities))
  
  
  if (!reference_team %in% team_names) {
    home_abilities[reference_team] <- 0
    away_abilities[reference_team] <- 0
    team_names <- c(team_names, reference_team)
  }
  
  
  home_ranks <- rank(-home_abilities)
  away_ranks <- rank(away_abilities)

  
  data.frame(Team = team_names, Home_Ability_Rank = home_ranks, Away_Ability_Rank = away_ranks)
}
extract_haf <- function(estimate, reference_team) {
  home_abilities <- -estimate[grepl("Home", names(estimate))]
  away_abilities <- -estimate[grepl("Away", names(estimate))]
  team_names <- gsub("Home|Away", "", names(home_abilities))
  
  
  if (!reference_team %in% team_names) {
    home_abilities[reference_team] <- 0
    away_abilities[reference_team] <- 0
    team_names <- c(team_names, reference_team)
  }
  home_advantage <- home_abilities - away_abilities
  

  haf_ranks <- rank(-home_advantage)
  
  data.frame(Team = team_names, HAF = home_advantage, HAF_ranks = haf_ranks)
}
pre_covid_haf <- extract_haf(pre_covid_train_fit$par, reference_team)
covid_haf <- extract_haf(covid_train_fit$par, reference_team)
post_covid_haf <- extract_haf(post_covid_train_fit$par, reference_team)

haf_table <- pre_covid_haf %>%
  rename('HAF (Pre-Covid)' = HAF, 'HAF Rank (Pre-Covid)' = HAF_ranks) %>%
  inner_join(covid_haf %>% rename('HAF (Covid)' = HAF, 'HAF Rank (Covid)' = HAF_ranks), by = "Team") %>%
  inner_join(post_covid_haf %>% rename('HAF (Post-Covid)' = HAF, 'HAF Rank (Post-Covid)' = HAF_ranks), by = "Team")
haf_table <- haf_table %>% arrange(`HAF Rank (Pre-Covid)`)

pre_covid_ranks <- extract_ranks(pre_covid_train_fit$par, reference_team)
covid_ranks <- extract_ranks(covid_train_fit$par, reference_team)
post_covid_ranks <- extract_ranks(post_covid_train_fit$par, reference_team)

final_rank_table <- pre_covid_ranks %>%
  rename('Home Ability Rank (Pre-Covid)' = Home_Ability_Rank, 'Away Ability Rank (Pre-Covid)' = Away_Ability_Rank) %>%
  inner_join(covid_ranks %>% rename('Home Ability Rank (Covid)' = Home_Ability_Rank, 'Away Ability Rank (Covid)' = Away_Ability_Rank), by = "Team") %>%
  inner_join(post_covid_ranks %>% rename('Home Ability Rank (Post-Covid)' = Home_Ability_Rank, 'Away Ability Rank (Post-Covid)' = Away_Ability_Rank), by = "Team")

#############################################################################
calculate_victory_probabilities <- function(data) {
  
  home_victories <- data %>% 
    filter(FTR == "H") %>% 
    count(Home) %>% 
    rename(Team = Home, Home_Victories = n)
  
  
  away_victories <- data %>% 
    filter(FTR == "A") %>% 
    count(Away) %>% 
    rename(Team = Away, Away_Victories = n)
  
  
  total_home_games <- data %>% count(Home) %>% rename(Team = Home, Total_Home_Games = n)
  total_away_games <- data %>% count(Away) %>% rename(Team = Away, Total_Away_Games = n)
  
  home_prob <- left_join(total_home_games, home_victories, by = "Team") %>%
    mutate(Home_Victory_Prob = ifelse(is.na(Home_Victories), 0, Home_Victories) / Total_Home_Games)
  away_prob <- left_join(total_away_games, away_victories, by = "Team") %>%
    mutate(Away_Victory_Prob = ifelse(is.na(Away_Victories), 0, Away_Victories) / Total_Away_Games)
  
  
  victory_prob <- left_join(home_prob, away_prob %>% select(Team, Away_Victory_Prob), by = "Team")
  
  return(victory_prob)
}
calculate_yearly_victory_probabilities <- function(data, periods) {
  yearly_victory_probs <- list()
  
  for (period in names(periods)) {
    years <- periods[[period]]
    yearly_victory_probs[[period]] <- lapply(years, function(year) {
      season_data <- filter(data, Season_End_Year == year)
      calculate_victory_probabilities(season_data)
    })
  }
  
  return(yearly_victory_probs)
}

calculate_period_average_victory_rank <- function(yearly_probs, team_appearances) {
  ranks <- list()
  
  for (period in names(yearly_probs)) {
    period_probs <- bind_rows(yearly_probs[[period]])
    
    period_probs <- period_probs %>%
      group_by(Team) %>%
      summarise(Home_Victory_Prob = mean(Home_Victory_Prob, na.rm = TRUE),
                Away_Victory_Prob = mean(Away_Victory_Prob, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(Team %in% team_appearances) %>%  
      mutate(Home_Victory_Rank = dense_rank(-Home_Victory_Prob),
             Away_Victory_Rank = dense_rank(-Away_Victory_Prob)) %>%
      arrange(Home_Victory_Rank)
    
    ranks[[period]] <- period_probs
  }
  
  return(ranks)
}
periods <- list(
  "Pre-Covid" = 2014:2018,
  "Covid" = 2020:2021,
  "Post-Covid" = 2023
)
team_appearances <- soccer %>%
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2023) %>%
  group_by(Home) %>%
  summarise(Distinct_Seasons = n_distinct(Season_End_Year)) %>%
  filter(Distinct_Seasons >= 10) %>%
  pull(Home)
yearly_victory_probs <- calculate_yearly_victory_probabilities(soccer, periods)
average_victory_ranks <- calculate_period_average_victory_rank(yearly_victory_probs,team_appearances)

final_rank_table <- final_rank_table %>%
  left_join(average_victory_ranks[["Pre-Covid"]] %>% 
              select(Team, Home_Victory_Rank, Away_Victory_Rank), by = "Team") %>%
  left_join(average_victory_ranks[["Covid"]] %>% 
              select(Team, Home_Victory_Rank, Away_Victory_Rank), by = "Team") %>%
  left_join(average_victory_ranks[["Post-Covid"]] %>% 
              select(Team, Home_Victory_Rank, Away_Victory_Rank), by = "Team") %>%
  rename(`Actual Home Victory Mean Rank (Pre-Covid)` = Home_Victory_Rank.x,
         `Actual Away Victory Mean Rank (Pre-Covid)` = Away_Victory_Rank.x,
         `Actual Home Victory Mean Rank (Covid)` = Home_Victory_Rank.y,
         `Actual Away Victory Mean Rank (Covid)` = Away_Victory_Rank.y,
         `Actual Home Victory Mean Rank (Post-Covid)` = Home_Victory_Rank,
         `Actual Away Victory Mean Rank (Post-Covid)` = Away_Victory_Rank)

final_rank_table <- final_rank_table %>%
  select(Team,
         "Home Ability Rank (Pre-Covid)", "Actual Home Victory Mean Rank (Pre-Covid)",
         "Away Ability Rank (Pre-Covid)", "Actual Away Victory Mean Rank (Pre-Covid)",
         "Home Ability Rank (Covid)", "Actual Home Victory Mean Rank (Covid)",
         "Away Ability Rank (Covid)", "Actual Away Victory Mean Rank (Covid)",
         "Home Ability Rank (Post-Covid)", "Actual Home Victory Mean Rank (Post-Covid)",
         "Away Ability Rank (Post-Covid)", "Actual Away Victory Mean Rank (Post-Covid)"
  )

periods_test <- list(
  "Pre-Covid" = 2019,
  "Covid" = 2022,
  "Post-Covid" = 2023
)
yearly_victory_test_probs <- calculate_yearly_victory_probabilities(soccer, periods_test)
prediction_victory_ranks <- calculate_period_average_victory_rank(yearly_victory_probs,team_appearances)
prediction_table <- final_rank_table %>%
  select(Team, "Home Ability Rank (Pre-Covid)","Away Ability Rank (Pre-Covid)", "Home Ability Rank (Covid)","Away Ability Rank (Covid)", "Home Ability Rank (Post-Covid)","Away Ability Rank (Post-Covid)") %>%
  left_join(prediction_victory_ranks[["Pre-Covid"]] %>% 
              select(Team, Home_Victory_Rank, Away_Victory_Rank), by = "Team") %>%
  left_join(prediction_victory_ranks[["Covid"]] %>% 
              select(Team, Home_Victory_Rank, Away_Victory_Rank), by = "Team") %>%
  left_join(prediction_victory_ranks[["Post-Covid"]] %>% 
              select(Team, Home_Victory_Rank, Away_Victory_Rank), by = "Team") %>%
  rename(`Actual Home Victory Rank (Pre-Covid)` = Home_Victory_Rank.x,
         `Actual Away Victory Rank (Pre-Covid)` = Away_Victory_Rank.x,
         `Actual Home Victory Rank (Covid)` = Home_Victory_Rank.y,
         `Actual Away Victory Rank (Covid)` = Away_Victory_Rank.y,
         `Actual Home Victory Rank (Post-Covid)` = Home_Victory_Rank,
         `Actual Away Victory Rank (Post-Covid)` = Away_Victory_Rank)
prediction_table <- prediction_table %>%
  select(Team,
         "Home Ability Rank (Pre-Covid)", "Actual Home Victory Rank (Pre-Covid)",
         "Away Ability Rank (Pre-Covid)", "Actual Away Victory Rank (Pre-Covid)",
         "Home Ability Rank (Covid)", "Actual Home Victory Rank (Covid)",
         "Away Ability Rank (Covid)", "Actual Away Victory Rank (Covid)",
         "Home Ability Rank (Post-Covid)", "Actual Home Victory Rank (Post-Covid)",
         "Away Ability Rank (Post-Covid)", "Actual Away Victory Rank (Post-Covid)"
  )
final_rank_table <- final_rank_table %>% arrange(`Home Ability Rank (Pre-Covid)`)
prediction_table <- prediction_table %>% arrange(`Home Ability Rank (Pre-Covid)`)

library(flextable)
haf_rank_table <- flextable(haf_table)
train_final_rank <- flextable(final_rank_table)
test_final_rank <- flextable(prediction_table)
print(train_final_rank)
print(test_final_rank)

save_as_docx(haf_rank_table, path = "BT-model2-haf.docx")
save_as_docx(train_final_rank, path = "BT-model2-train.docx")
save_as_docx(test_final_rank, path = "BT-model2-test.docx")