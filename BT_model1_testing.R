library(dplyr)
library(flextable)
library(officer)
library(numDeriv)


soccer <- read.csv('premier-league-matches.csv')
soccer$FTR <- factor(soccer$FTR, levels = c("A", "D", "H"), ordered = TRUE)
soccer$Home <- factor(soccer$Home)
soccer$Away <- factor(soccer$Away)
reference_team <- "Liverpool"
soccer$HomeAdvantage <- 1
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


soccer_filtered <- soccer_filtered %>%
  filter(Season_End_Year >= 2014 & Season_End_Year <= 2023)

######################################################################
log_likelihood <- function(params, data, reference) {
  
  delta <- params[28]
  home_advantage_pre <- params[29]
  home_advantage_covid <- params[30]
  home_advantage_post <- params[31]
  
  
  
  
  log_like <- 0
  
  for (i in 1:nrow(data)) {
    match <- data[i, ]
    home_suffix <- ifelse(match$pre_covid == 1, "pre",
                          ifelse(match$covid == 1, "covid", "post"))
    home_param_name <- ifelse(as.character(match$Home) == reference, reference, paste(as.character(match$Home), home_suffix, sep="_"))
    home_ability <- params[home_param_name]
    
    away_suffix <- ifelse(match$pre_covid == 1, "pre",
                          ifelse(match$covid == 1, "covid", "post"))
    away_param_name <- ifelse(as.character(match$Away) == reference, reference, paste(as.character(match$Away), away_suffix, sep="_"))
    home_ability <- ifelse(home_param_name==reference, 0, params[home_param_name])
    away_ability <- ifelse(away_param_name==reference, 0, params[away_param_name])
    
    home_advantage_value_pre <- ifelse(match$HomeAdvantage == 1 & match$pre_covid == 1, home_advantage_pre, 0)
    home_advantage_value_covid <- ifelse(match$HomeAdvantage == 1 & match$covid == 1, home_advantage_covid, 0)
    home_advantage_value_post <- ifelse(match$HomeAdvantage == 1 & match$post_covid == 1, home_advantage_post, 0)
    
    logit_D <- delta + home_ability - away_ability + home_advantage_value_pre + home_advantage_value_covid + home_advantage_value_post
    logit_A <- -delta + home_ability -  away_ability + home_advantage_value_pre + home_advantage_value_covid + home_advantage_value_post
    
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

teams <- unique(c(as.character(soccer_filtered$Home), as.character(soccer_filtered$Away)))
teams <- teams[teams != reference_team]
param_names <- unlist(lapply(teams, function(team) {
  c(paste(team, "pre", sep="_"), paste(team, "covid", sep="_"), paste(team, "post", sep="_"))
}))
initial_params <- c(rep(c(1, -1, 1), each = 9), delta = 1, home_advantage_pre = 0.5, home_advantage_covid = 0,home_advantage_post = 1)
names(initial_params) <- c(param_names, "delta", "home_advantage_pre","home_advantage_covid","home_advantage_post")
model_fit <- nlminb(initial_params,objective = log_likelihood, data = soccer_filtered, reference = reference_team, control=list(trace=1,x.tol=0.0000000001))
library(tidyr)
params_df <- data.frame(parameter = names(model_fit$par), value = -model_fit$par)


teams_df <- params_df %>%
  filter(!grepl("delta|home_advantage", parameter)) %>%
  separate(parameter, into = c("Team", "Period"), sep = "_") %>%
  spread(key = "Period", value = "value") %>%
  bind_rows(data.frame(Team = "Liverpool", pre = 0, covid = 0, post = 0)) %>%
  arrange(Team) %>%
  mutate(rank_pre = rank(-pre),
         rank_covid = rank(-covid),
         rank_post = rank(-post))


print(teams_df)
params_est <- model_fit$par
hessian <- hessian(log_likelihood, params_est, data = soccer_filtered, reference = reference_team)
cov_matrix <- solve(hessian)
contrast_matrix <- rbind(
  c(rep(0, 27), 0, -1, 1, 0),  # home_advantage_covid - home_advantage_pre
  c(rep(0, 27), 0, -1, 0, 1),  # home_advantage_post - home_advantage_pre
  c(rep(0, 27), 0, 0, -1, 1) # home_advantage_post - home_advantage_covid
)

contrast_est <- contrast_matrix %*% params_est
contrast_se <- sqrt(diag(contrast_matrix %*% cov_matrix %*% t(contrast_matrix)))
z_values <- contrast_est / contrast_se
p_values <- 2 * (1 - pnorm(abs(z_values)))
results <- data.frame(
  Contrast = c("home_advantage_covid - home_advantage_pre", "home_advantage_post - home_advantage_pre","home_advantage_post - home_advantage_covid"),
  Estimate = contrast_est,
  StdError = contrast_se,
  ZValue = z_values,
  PValue = p_values
)

print(results)
