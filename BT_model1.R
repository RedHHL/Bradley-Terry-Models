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
soccer$HomeAdvantage <- 1
reference_team <- "Liverpool"

team_appearances <- soccer %>% 
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2023) %>% 
  group_by(Home) %>% 
  summarise(Distinct_Seasons = n_distinct(Season_End_Year)) %>% 
  filter(Distinct_Seasons >= 10)
soccer_filtered <- soccer %>% 
  filter(Home %in% team_appearances$Home & Away %in% team_appearances$Home)
soccer_filtered$Home <- relevel(soccer_filtered$Home, ref = reference_team)
soccer_filtered$Away <- relevel(soccer_filtered$Away, ref = reference_team)
models <- list()
estimates <- list()
teams <- unique(c(as.character(soccer_filtered$Home), as.character(soccer_filtered$Away)))
pairings <- combn(teams, 2, simplify = FALSE)



split_2014_2019 <- split_data_by_year(soccer_filtered, 2014, 2018, 2019, 2019)
split_2020_2022 <- split_data_by_year(soccer_filtered, 2020, 2021, 2022, 2022)

#############################################################################

log_likelihood <- function(params, data, reference) {

  delta <- params[length(teams) + 1]
  home_advantage <- params[length(teams) + 2]


  
  log_like <- 0
  
  for (i in 1:nrow(data)) {
    match <- data[i, ]
    

    home_ability <- ifelse(as.character(match$Home) == reference, 0, params[as.character(match$Home)])
    away_ability <- ifelse(as.character(match$Away) == reference, 0, params[as.character(match$Away)])
    home_advantage_value <- ifelse(match$HomeAdvantage == 1, home_advantage, 0)
    
    logit_D <- delta + home_ability - away_ability + home_advantage_value
    logit_A <- -delta + home_ability - away_ability + home_advantage_value
    
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
teams <- unique(c(as.character(split_2014_2019$train$Home), as.character(split_2014_2019$train$Away)))
teams <- teams[teams != reference_team]
initial_params <- c(rep(1,9), delta = 1, home_advantage = 0.5)
names(initial_params) <- c(teams, "delta", "home_advantage")
pre_covid_train_fit <- nlminb(initial_params,objective = log_likelihood, data = split_2014_2019$train, reference = reference_team, control=list(trace=1,x.tol=0.0000000001))
covid_train_fit <- nlminb(initial_params,objective = log_likelihood, data = split_2020_2022$train, reference = reference_team, control=list(trace=1,x.tol=0.0000000001))
post_covid_train_fit <- nlminb(initial_params,objective = log_likelihood, data = filter(soccer_filtered, Season_End_Year == 2023), reference = reference_team, control=list(trace=1,x.tol=0.0000000001))

pre_covid_train_fit$par
covid_train_fit$par
post_covid_train_fit$par


#############################################################################
calculate_probabilities <- function(data, params, reference) {
  predictions <- data.frame(A = numeric(nrow(data)), D = numeric(nrow(data)), H = numeric(nrow(data)))
  delta <- params["delta"]
  home_advantage <- params["home_advantage"]
  for (i in 1:nrow(data)) {
    match <- data[i, ]
    home_ability <- ifelse(as.character(match$Home) == reference, 0, params[as.character(match$Home)])
    away_ability <- ifelse(as.character(match$Away) == reference, 0, params[as.character(match$Away)])
    home_advantage_value <- ifelse(match$HomeAdvantage == 1, home_advantage, 0)
    
  
    logit_D <- delta + home_ability - away_ability + home_advantage_value
    logit_A <- -delta + home_ability - away_ability + home_advantage_value
    
    cum_prob_D <- exp(logit_D) / (1 + exp(logit_D))
    cum_prob_A <- exp(logit_A) / (1 + exp(logit_A))
    prob_H <- 1 - cum_prob_D
    prob_D <- cum_prob_D - cum_prob_A
    prob_A <- cum_prob_A
    
    
    predictions[i, ] <- c(A = prob_A, D = prob_D, H = prob_H)
  }
  return(predictions)
}


calculate_scores <- function(data, params, reference, log_likelihood) {
  predicted_probs <- calculate_probabilities(data, params, reference)
  # Brier Score
  outcomes <- c("A", "D", "H")
  actual_probs <- matrix(as.numeric(unlist(lapply(data$FTR, function(x) outcomes == x))), nrow = nrow(data), byrow = TRUE)
  brier_score <- mean(rowSums((predicted_probs - actual_probs) ^ 2) / ncol(actual_probs))
  
  # Ranked Probability Score (RPS)
  predicted_cum_probs <- t(apply(predicted_probs, 1, cumsum))
  actual_cum_probs <- t(apply(actual_probs, 1, cumsum))
  rps <- mean(rowSums((predicted_cum_probs - actual_cum_probs) ^ 2) / (ncol(actual_probs) - 1))
  
  # AIC
  ll <- log_likelihood(params, data, reference)
  k <- length(params)
  aic <- 2 * k + 2 * ll
  bic <- log(length(data)) * k + 2 * ll
  
  return(list(BrierScore = brier_score, RankedProbabilityScore = rps, AIC = aic, BIC = bic, log_likelihood = ll))
}
pre_covid_scores <- calculate_scores(
  data = split_2014_2019$test,
  params = pre_covid_train_fit$par,
  reference = reference_team,
  log_likelihood = log_likelihood
)
covid_scores <- calculate_scores(
  data = split_2020_2022$test,
  params = covid_train_fit$par,
  reference = reference_team,
  log_likelihood = log_likelihood
)
pre_covid_scores 
covid_scores
#############################################################################
extract_and_rank <- function(estimate, reference_team) {
  abilities <- estimate[names(estimate)[1:9]]
  team_names <- names(abilities)
  abilities <- round(abilities, 4)
  if (!reference_team %in% team_names) {
    abilities[reference_team] <- 0
    team_names <- c(team_names, reference_team)
  }
  
  data.frame(Team = team_names, Ability_Parameter = -abilities, Ability_Rank = rank(abilities))
}
  
pre_covid_data <- extract_and_rank(pre_covid_train_fit$par,reference_team)
covid_data <- extract_and_rank(covid_train_fit$par,reference_team)
post_covid_data <- extract_and_rank(post_covid_train_fit$par,reference_team)

final_table <- pre_covid_data %>%
  rename('Ability Parameter (Pre-Covid)' = Ability_Parameter, 'Ability Rank (Pre-Covid)' = Ability_Rank) %>%
  inner_join(covid_data %>% rename('Ability Parameter (Covid)' = Ability_Parameter, 'Ability Rank (Covid)' = Ability_Rank), by = "Team") %>%
  inner_join(post_covid_data %>% rename('Ability Parameter (Post-Covid)' = Ability_Parameter, 'Ability Rank (Post-Covid)' = Ability_Rank), by = "Team")

#############################################################################

calculate_total_points <- function(data, win_points) {
  
  data <- data %>%
    mutate(Home_Points = 0, Away_Points = 0)
  
  data <- data %>%
    mutate(
      Home_Points = case_when(
        FTR == "H" ~ win_points,
        FTR == "D" ~ 1, 
        TRUE ~ 0 
      ),
      Away_Points = case_when(
        FTR == "A" ~ win_points, 
        FTR == "D" ~ 1,
        TRUE ~ 0  
      )
    )

  home_points <- data %>%
    group_by(Home) %>%
    summarize(Home_Total = sum(Home_Points)) %>%
    ungroup()

  
  away_points <- data %>%
    group_by(Away) %>%
    summarize(Away_Total = sum(Away_Points)) %>%
    ungroup()

  total_points <- home_points %>%
    rename(Team = Home) %>%
    full_join(away_points %>% rename(Team = Away), by = "Team") %>%
    mutate(
      Home_Total = coalesce(Home_Total, 0),
      Away_Total = coalesce(Away_Total, 0),
      Total_Points = Home_Total + Away_Total
    ) %>%
    arrange(desc(Total_Points))

  

  return(total_points)
  
}


calculate_seasonal_rankings <- function(season_data, win_points) {
  
  season_scores <- calculate_total_points(season_data, win_points)
  
  
  season_scores <- season_scores %>%
    arrange(desc(Total_Points)) %>%
    mutate(Rank = row_number()) %>%
    rename(Season_Rank = Rank)

  season_scores

}

calculate_rankings <- function(data, seasons, team_appearances, win_points) {
  
  seasonal_rankings_list <- lapply(seasons, function(season) {
    season_data <- filter(data, Season_End_Year == season)
    calculate_seasonal_rankings(season_data, win_points)
  })
  
  seasonal_rankings <- bind_rows(seasonal_rankings_list)

  
  
  average_rankings <- seasonal_rankings %>%
    group_by(Team) %>%
    summarise(Average_Rank = mean(Season_Rank, na.rm = TRUE)) %>%
    ungroup()
  
  final_rankings <- average_rankings %>%
    filter(Team %in% team_appearances) %>%
    arrange(Average_Rank) %>%
    mutate(Final_Rank = row_number())

  
  final_rankings
}


pre_covid_years <- 2014:2018
covid_years <- 2020:2021
post_covid_year <- 2023


team_appearances <- soccer %>%
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2023) %>%
  group_by(Home) %>%
  summarise(Distinct_Seasons = n_distinct(Season_End_Year)) %>%
  filter(Distinct_Seasons >= 10) %>%
  pull(Home)

rankings_epl_pre_covid <- calculate_rankings(soccer, pre_covid_years, team_appearances, 3)
rankings_epl_covid <- calculate_rankings(soccer, covid_years, team_appearances, 3)
rankings_epl_post_covid <- calculate_rankings(soccer, list(post_covid_year), team_appearances, 3)


rankings_model_pre_covid <- calculate_rankings(soccer, pre_covid_years, team_appearances, 2)
rankings_model_covid <- calculate_rankings(soccer, covid_years, team_appearances, 2)
rankings_model_post_covid <- calculate_rankings(soccer, list(post_covid_year), team_appearances, 2)



rankings_model_pre_covid_filtered <- calculate_rankings(soccer_filtered, pre_covid_years, team_appearances, 2)
rankings_model_covid_filtered <- calculate_rankings(soccer_filtered, covid_years, team_appearances, 2)
rankings_model_post_covid_filtered <- calculate_rankings(soccer_filtered, list(post_covid_year), team_appearances, 2)




final_table <- final_table %>%
  left_join(rankings_epl_pre_covid %>% dplyr::select(Team, Actual_Mean_Rank_Pre_Covid = Final_Rank), by = "Team") %>%
  left_join(rankings_epl_covid %>% dplyr::select(Team, Actual_Mean_Rank_Covid = Final_Rank), by = "Team") %>%
  left_join(rankings_epl_post_covid %>% dplyr::select(Team, Actual_Mean_Rank_Post_Covid = Final_Rank), by = "Team") %>%
  left_join(rankings_model_pre_covid %>% dplyr::select(Team, Actual_Mean_Rank_Model_Pre_Covid = Final_Rank), by = "Team") %>%
  left_join(rankings_model_covid %>% dplyr::select(Team, Actual_Mean_Rank_Model_Covid = Final_Rank), by = "Team") %>%
  left_join(rankings_model_post_covid %>% dplyr::select(Team, Actual_Mean_Rank_Model_Post_Covid = Final_Rank), by = "Team") %>%
  left_join(rankings_model_pre_covid_filtered %>% dplyr::select(Team, Actual_Mean_Rank_Model_Pre_Covid_Filtered = Final_Rank), by = "Team") %>%
  left_join(rankings_model_covid_filtered %>% dplyr::select(Team, Actual_Mean_Rank_Model_Covid_Filtered = Final_Rank), by = "Team") %>%
  left_join(rankings_model_post_covid_filtered %>% dplyr::select(Team, Actual_Mean_Rank_Model_Post_Covid_Filtered = Final_Rank), by = "Team")



final_table <- final_table %>%
  dplyr::select(Team, 
         "Ability Parameter (Pre-Covid)", "Ability Rank (Pre-Covid)", Actual_Mean_Rank_Pre_Covid, Actual_Mean_Rank_Model_Pre_Covid, Actual_Mean_Rank_Model_Pre_Covid_Filtered,
         "Ability Parameter (Covid)", "Ability Rank (Covid)", Actual_Mean_Rank_Covid, Actual_Mean_Rank_Model_Covid, Actual_Mean_Rank_Model_Covid_Filtered,
         "Ability Parameter (Post-Covid)", "Ability Rank (Post-Covid)", Actual_Mean_Rank_Post_Covid, Actual_Mean_Rank_Model_Post_Covid, Actual_Mean_Rank_Model_Post_Covid_Filtered)


names(final_table) <- c("Team", "Ability Parameter (Pre-Covid)", "Ability Rank (Pre-Covid)",
                        "Actual Mean Rank in EPL (Pre-Covid)", "Mean Rank Using 012 scoring system(Pre-Covid)","Mean Rank Using 012 scoring system and filtered data(Pre-Covid)", "Ability Parameter (Covid)", 
                        "Ability Rank (Covid)", "Actual Mean Rank in EPL(Covid)", "Mean Rank Using 012 scoring system(Covid)","Mean Rank Using 012 scoring system and filtered data(Covid)",
                        "Ability Parameter (Post-Covid)", "Ability Rank (Post-Covid)", 
                        "Actual Mean Rank in EPL(Post-Covid)","Mean Rank Using 012 scoring system(Post-Covid)", "Mean Rank Using 012 scoring system and filtered data(Post-Covid)")
final_table <- final_table %>% arrange(`Ability Rank (Pre-Covid)`)
pre_covid_train_table <- final_table %>% 
  dplyr::select(Team, 
         "Ability Parameter (Pre-Covid)", 
         "Ability Rank (Pre-Covid)",
         "Actual Mean Rank in EPL (Pre-Covid)", 
         "Mean Rank Using 012 scoring system(Pre-Covid)",
         "Mean Rank Using 012 scoring system and filtered data(Pre-Covid)")

covid_train_table <- final_table %>% 
  dplyr::select(Team, 
         "Ability Parameter (Covid)", 
         "Ability Rank (Covid)", 
         "Actual Mean Rank in EPL(Covid)", 
         "Mean Rank Using 012 scoring system(Covid)",
         "Mean Rank Using 012 scoring system and filtered data(Covid)")

post_covid_train_table <- final_table %>% 
  dplyr::select(Team, 
         "Ability Parameter (Post-Covid)", 
         "Ability Rank (Post-Covid)", 
         "Actual Mean Rank in EPL(Post-Covid)",
         "Mean Rank Using 012 scoring system(Post-Covid)", 
         "Mean Rank Using 012 scoring system and filtered data(Post-Covid)")





rankings_test_epl_pre_covid <- calculate_rankings(soccer, 2019, team_appearances, 3)
rankings_test_epl_covid <- calculate_rankings(soccer, 2021, team_appearances, 3)
rankings_test_epl_post_covid <- calculate_rankings(soccer, list(post_covid_year), team_appearances, 3)


rankings_test_model_pre_covid <- calculate_rankings(soccer, 2019, team_appearances, 2)
rankings_test_model_covid <- calculate_rankings(soccer, 2021, team_appearances, 2)
rankings_test_model_post_covid <- calculate_rankings(soccer, list(post_covid_year), team_appearances, 2)

rankings_test_model_pre_covid_filtered <- calculate_rankings(soccer_filtered, 2019, team_appearances, 2)
rankings_test_model_covid_filtered <- calculate_rankings(soccer_filtered, 2021, team_appearances, 2)
rankings_test_model_post_covid_filtered <- calculate_rankings(soccer_filtered, list(post_covid_year), team_appearances, 2)

prediction_table <- final_table %>%
  dplyr::select(Team, "Ability Rank (Pre-Covid)", "Ability Rank (Covid)", "Ability Rank (Post-Covid)") %>%
  left_join(rankings_test_epl_pre_covid %>% rename("Actual Rank in EPL (Pre-Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_model_pre_covid %>% rename("Rank Using 012 scoring system(Pre-Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_model_pre_covid_filtered %>% rename("Rank Using 012 scoring system and filtered data(Pre-Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_epl_covid %>% rename("Actual Rank in EPL (Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_model_covid %>% rename("Rank Using 012 scoring system(Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_model_covid_filtered %>% rename("Rank Using 012 scoring system and filtered data(Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_epl_post_covid %>% rename("Actual Rank in EPL (Post-Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_model_post_covid %>% rename("Rank Using 012 scoring system(Post-Covid)" = Final_Rank), by = "Team") %>%
  left_join(rankings_test_model_post_covid_filtered %>% rename("Rank Using 012 scoring system and filtered data(Post-Covid)" = Final_Rank), by = "Team") %>%
  arrange(`Ability Rank (Pre-Covid)`) 
prediction_table <- prediction_table %>% arrange(`Ability Rank (Pre-Covid)`)

prediction_table <- prediction_table %>%
  dplyr::select(
    "Team", 
    "Ability Rank (Pre-Covid)", "Actual Rank in EPL (Pre-Covid)", "Rank Using 012 scoring system(Pre-Covid)", "Rank Using 012 scoring system and filtered data(Pre-Covid)",
    "Ability Rank (Covid)", "Actual Rank in EPL (Covid)", "Rank Using 012 scoring system(Covid)", "Rank Using 012 scoring system and filtered data(Covid)",
    "Ability Rank (Post-Covid)", "Actual Rank in EPL (Post-Covid)", "Rank Using 012 scoring system(Post-Covid)","Rank Using 012 scoring system and filtered data(Post-Covid)"
  )

pre_covid_test_table <- prediction_table %>% 
  dplyr::select("Team", 
         "Ability Rank (Pre-Covid)",
         "Actual Rank in EPL (Pre-Covid)", 
         "Rank Using 012 scoring system(Pre-Covid)",
         "Rank Using 012 scoring system and filtered data(Pre-Covid)")

covid_test_table <- prediction_table %>% 
  dplyr::select("Team",
         "Ability Rank (Covid)", 
         "Actual Rank in EPL (Covid)", 
         "Rank Using 012 scoring system(Covid)",
         "Rank Using 012 scoring system and filtered data(Covid)")

post_covid_test_table <- prediction_table %>% 
  dplyr::select("Team",
         "Ability Rank (Post-Covid)", 
         "Actual Rank in EPL (Post-Covid)",
         "Rank Using 012 scoring system(Post-Covid)", 
         "Rank Using 012 scoring system and filtered data(Post-Covid)")




flextable_pre_covid_train <- flextable(pre_covid_train_table)
flextable_covid_train <- flextable(covid_train_table)
flextable_post_covid_train <- flextable(post_covid_train_table)
flextable_output <- flextable(final_table)
flextable_final <- flextable(final_table)
flextable_prediction <- flextable(prediction_table)
flextable_pre_covid_test <- flextable(pre_covid_test_table)
flextable_covid_test <- flextable(covid_test_table)
flextable_post_covid_test <- flextable(post_covid_test_table)


save_as_docx(flextable_pre_covid_train, path = "BT-model1-train-Pre-Covid.docx")
save_as_docx(flextable_covid_train, path = "BT-model1-train-Covid.docx")
save_as_docx(flextable_post_covid_train, path = "BT-model1-train-Post-Covid.docx")
save_as_docx(flextable_pre_covid_test, path = "BT-model1-test-Pre-Covid.docx")
save_as_docx(flextable_covid_test, path = "BT-model1-test-Covid.docx")
save_as_docx(flextable_post_covid_test, path = "BT-model1-test-Post-Covid.docx")
print(flextable_pre_covid_train)
print(flextable_covid_train)
print(flextable_post_covid_train)
print(flextable_pre_covid_test)
print(flextable_covid_test)
print(flextable_post_covid_test)