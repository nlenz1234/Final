thresholds <- list()

for (i in 1:100) {
  set.seed(888 + i)
  trainIndex <- createDataPartition(
    data_model$occupancy_classifier,
    p = 0.7,
    list = FALSE,
    times = 1
  )
  
  dataTrain <- data_model[trainIndex, ]
  dataTest <- data_model[-trainIndex, ]
  
  from_sum <- dataTrain %>%
    group_by(from) %>%
    summarize(from_count = n(),
              from_sum = sum(occupancy_classifier, na.rm = TRUE))
  
  to_sum <- dataTrain %>%
    group_by(to) %>%
    summarize(to_count = n(),
              to_sum = sum(occupancy_classifier, na.rm = TRUE))
  
  stations_train <- stations_sf[, c(1, 12)] %>%
    left_join(from_sum, by = c("code" = "from")) %>%
    left_join(to_sum, by = c("code" = "to")) %>%
    replace(is.na(.), 0) %>%
    mutate(overall_count = from_count + to_count,
           overall_sum = from_sum + to_sum)
  
  coords <- st_coordinates(stations_train)
  
  knn_nb <- knn2nb(knearneigh(coords, k = 5))
  
  neighbor_results <- lapply(1:length(knn_nb), function(i) {
    neighbors <- knn_nb[[i]]
    
    overall_count_sum <- sum(stations_train$overall_count[neighbors], na.rm = TRUE)
    overall_sum_sum <- sum(stations_train$overall_sum[neighbors], na.rm = TRUE)
    
    return(
      data.frame(
        overall_count_sum = overall_count_sum,
        overall_sum_sum = overall_sum_sum
      )
    )
  })
  
  neighbor_results_df <- do.call(rbind, neighbor_results)
  
  stations_train <- stations_train %>%
    bind_cols(neighbor_results_df) %>%
    mutate(neighbor = overall_sum_sum / overall_count_sum) %>%
    st_drop_geometry()
  
  dataTrain <- dataTrain %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(from_neighbor = neighbor),
      by = c("from" = "code")
    ) %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(to_neighbor = neighbor),
      by = c("to" = "code")
    ) %>%
    mutate(OD_neighbor = (from_neighbor + to_neighbor) / 2)
  
  dataTrain <- dataTrain %>%
    mutate(
      mapped_date = as.Date(start_date) + ((as.numeric(
        as.Date(datetime_combined) - as.Date(start_date)
      )) %% 7),
      mapped_datetime = as.POSIXct(paste(
        mapped_date, format(datetime_combined, "%H:%M:%S")
      ), format = "%Y-%m-%d %H:%M:%S")
    )
  
  datetime <- dataTrain %>%
    dplyr::select(occupancy_classifier, datetime_combined, mapped_datetime) %>%
    arrange(mapped_datetime) %>%
    slice(1:20) %>%
    mutate(mapped_datetime = mapped_datetime + days(7)) %>%
    bind_rows(
      dataTrain %>%
        dplyr::select(occupancy_classifier, datetime_combined, mapped_datetime) %>%
        arrange(mapped_datetime) %>%
        slice((n() - 19):n()) %>%
        mutate(mapped_datetime = mapped_datetime - days(7))
    ) %>%
    bind_rows(
      dataTrain %>%
        dplyr::select(occupancy_classifier, datetime_combined, mapped_datetime)
    )
  
  datetime <- datetime %>%
    arrange(mapped_datetime)
  
  dataTrain <- dataTrain %>%
    mutate(Adjacant = map_dbl(1:n(), function(i) {
      distances <- abs(as.numeric(datetime$mapped_datetime - mapped_datetime[i]))
      nearest_indices <- order(distances)[1:20]
      mean(datetime$occupancy_classifier[nearest_indices], na.rm = TRUE)
    }))
  
  model <- glm(
    occupancy_classifier ~ type + Temperature + OD_passenger + Frequency + OD_neighbor + Adjacant,
    data = dataTrain,
    family = "binomial" (link = "logit")
  )
  
  dataTest <- dataTest %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(from_neighbor = neighbor),
      by = c("from" = "code")
    ) %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(to_neighbor = neighbor),
      by = c("to" = "code")
    ) %>%
    mutate(OD_neighbor = (from_neighbor + to_neighbor) / 2)
  
  dataTest <- dataTest %>%
    mutate(
      mapped_date = as.Date(start_date) + ((as.numeric(
        as.Date(datetime_combined) - as.Date(start_date)
      )) %% 7),
      mapped_datetime = as.POSIXct(paste(
        mapped_date, format(datetime_combined, "%H:%M:%S")
      ), format = "%Y-%m-%d %H:%M:%S")
    )
  
  dataTest <- dataTest %>%
    mutate(Adjacant = map_dbl(1:n(), function(i) {
      distances <- abs(as.numeric(datetime$mapped_datetime - mapped_datetime[i]))
      nearest_indices <- order(distances)[1:20]
      mean(datetime$occupancy_classifier[nearest_indices], na.rm = TRUE)
    }))
  
  dataTest <- dataTest %>%
    mutate(occupancy_classifier = factor(
      occupancy_classifier,
      levels = c(0, 1),
      labels = c("Low", "High")
    ))
  
  testProbs <- data.frame(
    Outcome = as.factor(dataTest$occupancy_classifier),
    Probs = predict(model, dataTest, type = "response")
  )
  
  iterateThresholds <- function(data) {
    x = .01
    all_prediction <- data.frame()
    while (x <= 1) {
      this_prediction <-
        testProbs %>%
        mutate(predOutcome = ifelse(Probs > x, "High", "Low")) %>%
        count(predOutcome, Outcome) %>%
        summarize(
          True_Negative = sum(n[predOutcome == "Low" & Outcome == "Low"]),
          True_Positive = sum(n[predOutcome == "High" &
                                  Outcome == "High"]),
          False_Negative = sum(n[predOutcome == "Low" &
                                   Outcome == "High"]),
          False_Positive = sum(n[predOutcome == "High" &
                                   Outcome == "Low"])
        ) %>%
        gather(Variable, Count) %>%
        mutate(
          Revenue =
            ifelse(
              Variable == "True_Negative",
              0 * Count,
              ifelse(
                Variable == "True_Positive",
                2 * Count,
                ifelse(
                  Variable == "False_Negative",
                  0 * Count,
                  ifelse(Variable == "False_Positive", (-1) * Count, 0)
                )
              )
            ),
          Threshold = x
        )
      
      all_prediction <- rbind(all_prediction, this_prediction)
      x <- x + .01
    }
    return(all_prediction)
  }
  
  whichThreshold <- iterateThresholds(testProbs)
  
  whichThreshold_revenue <-
    whichThreshold %>%
    group_by(Threshold) %>%
    summarize(Revenue = sum(Revenue))
  
  optimal_threshold <- whichThreshold_revenue %>%
    filter(Revenue == max(Revenue)) %>%
    pull(Threshold)
  
  max_revenue <- max(whichThreshold_revenue$Revenue)
  
  thresholds[[i]] <- list(max_revenue = max_revenue)
  
  # cat("Iteration:", i, "Optimal Threshold:", optimal_threshold, "\n")
}

avg_optimal <- mean(sapply(thresholds, function(res) res$max_revenue[1]), na.rm = T)
# cat("Average Optimal Threshold over 100 runs:", avg_optimal_threshold, "\n")





thresholds <- list()

for (i in 1:100) {
  set.seed(888 + i)
  trainIndex <- createDataPartition(
    data_model$occupancy_classifier,
    p = 0.7,
    list = FALSE,
    times = 1
  )
  
  
  rare_levels <- data_model %>%
    filter(hour == "3" | hour == "4" | hour == "2")
  remaining_data <- data_model %>%
    filter(!(hour == "3" | hour == "4" | hour == "2"))
  
  trainIndex <- createDataPartition(data_model$occupancy_classifier, p = 0.7, list = FALSE, times = 1)
  
  dataTrain <- bind_rows(rare_levels, remaining_data[trainIndex, ])
  dataTest <- remaining_data[-trainIndex, ]
  
  
  
  from_sum <- dataTrain %>%
    group_by(from) %>%
    summarize(from_count = n(),
              from_sum = sum(occupancy_classifier, na.rm = TRUE))
  
  to_sum <- dataTrain %>%
    group_by(to) %>%
    summarize(to_count = n(),
              to_sum = sum(occupancy_classifier, na.rm = TRUE))
  
  stations_train <- stations_sf[, c(1, 12)] %>%
    left_join(from_sum, by = c("code" = "from")) %>%
    left_join(to_sum, by = c("code" = "to")) %>%
    replace(is.na(.), 0) %>%
    mutate(overall_count = from_count + to_count,
           overall_sum = from_sum + to_sum)
  
  coords <- st_coordinates(stations_train)
  
  knn_nb <- knn2nb(knearneigh(coords, k = 5))
  
  neighbor_results <- lapply(1:length(knn_nb), function(i) {
    neighbors <- knn_nb[[i]]
    
    overall_count_sum <- sum(stations_train$overall_count[neighbors], na.rm = TRUE)
    overall_sum_sum <- sum(stations_train$overall_sum[neighbors], na.rm = TRUE)
    
    return(
      data.frame(
        overall_count_sum = overall_count_sum,
        overall_sum_sum = overall_sum_sum
      )
    )
  })
  
  neighbor_results_df <- do.call(rbind, neighbor_results)
  
  stations_train <- stations_train %>%
    bind_cols(neighbor_results_df) %>%
    mutate(neighbor = overall_sum_sum / overall_count_sum) %>%
    st_drop_geometry()
  
  dataTrain <- dataTrain %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(from_neighbor = neighbor),
      by = c("from" = "code")
    ) %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(to_neighbor = neighbor),
      by = c("to" = "code")
    ) %>%
    mutate(OD_neighbor = (from_neighbor + to_neighbor) / 2)
  
  dataTrain <- dataTrain %>%
    mutate(
      mapped_date = as.Date(start_date) + ((as.numeric(
        as.Date(datetime_combined) - as.Date(start_date)
      )) %% 7),
      mapped_datetime = as.POSIXct(paste(
        mapped_date, format(datetime_combined, "%H:%M:%S")
      ), format = "%Y-%m-%d %H:%M:%S")
    )
  
  datetime <- dataTrain %>%
    dplyr::select(occupancy_classifier, datetime_combined, mapped_datetime) %>%
    arrange(mapped_datetime) %>%
    slice(1:20) %>%
    mutate(mapped_datetime = mapped_datetime + days(7)) %>%
    bind_rows(
      dataTrain %>%
        dplyr::select(occupancy_classifier, datetime_combined, mapped_datetime) %>%
        arrange(mapped_datetime) %>%
        slice((n() - 19):n()) %>%
        mutate(mapped_datetime = mapped_datetime - days(7))
    ) %>%
    bind_rows(
      dataTrain %>%
        dplyr::select(occupancy_classifier, datetime_combined, mapped_datetime)
    )
  
  datetime <- datetime %>%
    arrange(mapped_datetime)
  
  dataTrain <- dataTrain %>%
    mutate(Adjacant = map_dbl(1:n(), function(i) {
      distances <- abs(as.numeric(datetime$mapped_datetime - mapped_datetime[i]))
      nearest_indices <- order(distances)[1:20]
      mean(datetime$occupancy_classifier[nearest_indices], na.rm = TRUE)
    }))
  
  model <- glm(
    occupancy_classifier ~ type + hour + dotw,
    data = dataTrain,
    family = "binomial" (link = "logit")
  )
  
  dataTest <- dataTest %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(from_neighbor = neighbor),
      by = c("from" = "code")
    ) %>%
    left_join(
      stations_train %>% select(code, neighbor) %>% rename(to_neighbor = neighbor),
      by = c("to" = "code")
    ) %>%
    mutate(OD_neighbor = (from_neighbor + to_neighbor) / 2)
  
  dataTest <- dataTest %>%
    mutate(
      mapped_date = as.Date(start_date) + ((as.numeric(
        as.Date(datetime_combined) - as.Date(start_date)
      )) %% 7),
      mapped_datetime = as.POSIXct(paste(
        mapped_date, format(datetime_combined, "%H:%M:%S")
      ), format = "%Y-%m-%d %H:%M:%S")
    )
  
  dataTest <- dataTest %>%
    mutate(Adjacant = map_dbl(1:n(), function(i) {
      distances <- abs(as.numeric(datetime$mapped_datetime - mapped_datetime[i]))
      nearest_indices <- order(distances)[1:20]
      mean(datetime$occupancy_classifier[nearest_indices], na.rm = TRUE)
    }))
  
  dataTest <- dataTest %>%
    mutate(occupancy_classifier = factor(
      occupancy_classifier,
      levels = c(0, 1),
      labels = c("Low", "High")
    ))
  
  testProbs <- data.frame(
    Outcome = as.factor(dataTest$occupancy_classifier),
    Probs = predict(model, dataTest, type = "response")
  )
  
  iterateThresholds <- function(data) {
    x = .01
    all_prediction <- data.frame()
    while (x <= 1) {
      this_prediction <-
        testProbs %>%
        mutate(predOutcome = ifelse(Probs > x, "High", "Low")) %>%
        count(predOutcome, Outcome) %>%
        summarize(
          True_Negative = sum(n[predOutcome == "Low" & Outcome == "Low"]),
          True_Positive = sum(n[predOutcome == "High" &
                                  Outcome == "High"]),
          False_Negative = sum(n[predOutcome == "Low" &
                                   Outcome == "High"]),
          False_Positive = sum(n[predOutcome == "High" &
                                   Outcome == "Low"])
        ) %>%
        gather(Variable, Count) %>%
        mutate(
          Revenue =
            ifelse(
              Variable == "True_Negative",
              0 * Count,
              ifelse(
                Variable == "True_Positive",
                2 * Count,
                ifelse(
                  Variable == "False_Negative",
                  0 * Count,
                  ifelse(Variable == "False_Positive", (-1) * Count, 0)
                )
              )
            ),
          Threshold = x
        )
      
      all_prediction <- rbind(all_prediction, this_prediction)
      x <- x + .01
    }
    return(all_prediction)
  }
  
  whichThreshold <- iterateThresholds(testProbs)
  
  whichThreshold_revenue <-
    whichThreshold %>%
    group_by(Threshold) %>%
    summarize(Revenue = sum(Revenue))
  
  optimal_threshold <- whichThreshold_revenue %>%
    filter(Revenue == max(Revenue)) %>%
    pull(Threshold)
  
  max_revenue <- max(whichThreshold_revenue$Revenue)
  
  thresholds[[i]] <- list(max_revenue = max_revenue)
  
  # cat("Iteration:", i, "Optimal Threshold:", optimal_threshold, "\n")
}

avg_optimal <- mean(sapply(thresholds, function(res) res$max_revenue[1]), na.rm = T)
# cat("Average Optimal Threshold over 100 runs:", avg_optimal_threshold, "\n")
results <- list()


