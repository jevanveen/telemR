#needed fixes/features
#summary study.
#in graph, last dark box not drawing
#wanted fixes/features

#' Combine factors in a tidy telem giving a new column
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi. MUST HAVE grouping meta data in Xover
#' @return a one day view of uncrossed telem data
#' @export
x_over_telem <- function(tidy_telem, day_1_start, day_2_start, xover_var = "Treated", xover_pattern_1 = "A", xover_pattern_2 = "B"){
  tryCatch({
    day1 <- trim_telem(tidy_telem = tidy_telem, start_time = day_1_start, end_time = (as.POSIXct(day_1_start) + lubridate::days(1))) %>%
      mutate(Xover = ifelse(grepl(xover_pattern_1, Xover), xover_var, "Control"))
    day2 <- trim_telem(tidy_telem = tidy_telem, start_time = day_2_start, end_time = (as.POSIXct(day_2_start) + lubridate::days(1))) %>%
      mutate(Xover = ifelse(grepl(xover_pattern_2, Xover), xover_var, "Control"))
    ret <- combine_telem(day1, day2)
    ret$Xover <- as.factor(ret$Xover)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(ret)
}
#' Combine factors in a tidy telem giving a new column
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param telem_var_1 first telem var to combine
#' @param telem_var_2 second telem var to combine
#' @return a  tidy telemetry tibble with a new column "x_factor"
#' @export
x_factor_telem <- function(tidy_telem, telem_var_1, telem_var_2){
  tryCatch({
    ret <- tidy_telem %>%
      mutate(x_factor = interaction(tidy_telem[[telem_var_1]], tidy_telem[[telem_var_2]]))
    ret$x_factor <- as.factor(ret$x_factor)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(ret)
}
#' Take two tidy telem files and align(combine) them, without collpasing.
#' Note that new dates for tidy_telem_2 will be inaccurate.
#' If animals in the two experiments have the same number and meta data,
#' their values will be averaged by default
#'
#' @param tidy_telem_1 a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param tidy_telem_2 a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param align_date_1 date in tidy_telem_1 to align with date in tidy_telem_2
#' @param align_date_2 date in tidy_telem_2 to align with date in tidy_telem_1
#' @return a collpased tidy telemetry tibble
#' @export
align_telem <- function(tidy_telem_1, tidy_telem_2, align_date_1 = NULL, align_date_2 = NULL){
  tryCatch({
    if(is.null(align_date_1)){
      align_date_1 <- first_time(tidy_telem_1$Time) %>% lubridate::date()
    }
    if(is.null(align_date_2)){
      align_date_2 <- first_time(tidy_telem_2$Time) %>% lubridate::date()
    }
      t1 <- as.POSIXct(paste0(align_date_1, " ", "23:55:55"))
      t2 <- as.POSIXct(paste0(align_date_2, " ", "23:55:55"))
      tdiff <- t2 - t1
      tidy_telem_2$Time <- tidy_telem_2$Time - tdiff
      ret <- bind_rows(tidy_telem_1, tidy_telem_2)
      if("Counts" %in% colnames(ret)){
        ret <- ret %>%
          group_by_at(setdiff(names(ret), c("DegC", "Counts"))) %>%
          summarise(Counts = mean(Counts, na.rm = T), DegC = mean(DegC, na.rm = T))
      } else {
        ret <- ret %>%
          group_by_at(setdiff(names(ret), c("DegC"))) %>%
          summarise(DegC = mean(DegC, na.rm = T))
      }

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  ret <- ungroup(ret)
  return(ret)
}
#' Take two tidy telem files and combine them,
#'
#' @param tidy_telem_1 a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param tidy_telem_2 a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @return a collpased tidy telemetry tibble
#' @export
combine_telem <- function(tidy_telem_1, tidy_telem_2){
  tryCatch({
    if(is.POSIXct(tidy_telem_1$Time)){
      c1 <- collapse_telem(tidy_telem_1)
    }
    if(is.POSIXct(tidy_telem_2$Time)){
      c2 <- collapse_telem(tidy_telem_2)
    }
    ret <- bind_rows(c1,c2)
    if("Counts" %in% colnames(ret)){
      ret <- ret %>%
        group_by_at(setdiff(names(ret), c("DegC", "Counts"))) %>%
        summarise(Counts = mean(Counts, na.rm = T), DegC = mean(DegC, na.rm = T))
    } else {
      ret <- ret %>%
        group_by_at(setdiff(names(ret), c("DegC"))) %>%
        summarise(DegC = mean(DegC, na.rm = T))
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  ret <- ungroup(ret)
  return(ret)
}
#' Read .asc files produced by starr telemetry and put them in tidy tibbles
#'
#' @param raw_asc .asc file produced by starr telemetry system. Use "new format" and "date-time in left column only" when exporting .asc
#' @param meta_data_group data frame or tibble containing sample group information
#' @param meta_data_xover data frame or tibble containing sample corssover groups
#' @param meta_data_sex data frame or tibble containing sample sexes
#' @param n_animals number of animals
#' @param n_measurements number of measurments taken
#' @param trim_na logical: should function remove time points containing NA values?
#' @param trim_bad_probe logical: should function remove bad probes?
#' @return a tidy telemetry tibble (tidy_telem)
#' @export
read_starr <- function(raw_asc, meta_data_group, meta_data_xover = NULL, meta_data_sex = NULL, trim_na = F, trim_bad_probe = F){
  tryCatch({
    skip <- grep(pattern = "YY/MM/DD", x = read.delim(file = raw_asc, blank.lines.skip = F, header = F)$V1)
    temp <- read_tsv(file = raw_asc, skip = (skip-1))
    if(trim_na == T){temp <- na.omit(temp)}
    temp[[1]] <- as.POSIXct(temp[[1]])
    colnames(temp)[1] <- "Time"
    suppressWarnings(
      samples <- read_tsv(file = raw_asc) %>%
        t() %>%
        grep(pattern = "Animal ID:", value = T) %>%
        sub(pattern = "Animal ID: ", replacement = "") %>%
        unique()
    )
    if(is.null(meta_data_sex)){meta_data_sex <- tibble("Unknown" = samples)}
    if(is.null(meta_data_xover)){meta_data_xover <- tibble("NA" = samples)}
    tlist <- list()
    for(i in samples){
      tlist[[i]] <- dplyr::select(temp, contains(i), contains("Time"))
      if(grepl("Cnts", colnames(tlist[[i]][,1]))){colnames(tlist[[i]])[1] <- "Counts"}
      if(grepl("Cnts", colnames(tlist[[i]][,2]))){colnames(tlist[[i]])[2] <- "Counts"}
      if(grepl("Deg. C", colnames(tlist[[i]][,1]))){colnames(tlist[[i]])[1] <- "DegC"}
      if(grepl("Deg. C", colnames(tlist[[i]][,2]))){colnames(tlist[[i]])[2] <- "DegC"}
      tlist[[i]][,"Mouse"] <- i
      tlist[[i]][,"Group"] <- colnames(meta_data_group[grep(i,meta_data_group)])
      tlist[[i]][,"Xover"] <- colnames(meta_data_xover[grep(i,meta_data_xover)])
      tlist[[i]][,"Sex"] <- colnames(meta_data_sex[grep(i,meta_data_sex)])
    }
    collapse <- bind_rows(tlist)
    collapse$Mouse <- as.factor(collapse$Mouse)
    collapse$Group <- as.factor(collapse$Group)
    collapse$Xover <- as.factor(collapse$Xover)
    collapse$Sex <- as.factor(collapse$Sex)
    collapse$Time <- collapse$Time %>%
      + lubridate::years(2000)
    if(trim_bad_probe == T){collpase <- collpase %>% trim_bad_probe()}

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(collapse)

}
#' Read a folder of xlsx files produced by star-oddi telemetry and put them together in a tidy tibble
#'
#' @param folder a folder containing xlsx files produced by star-oddi emitters
#' @param meta_data data frame or tibble containing sample group information
#' @param trim_bad_probe logical: should function remove bad probes? Particularly useful for probes attached to mouse, rather than implanted
#' @return a tidy telemetry tibble (tidy_telem)
#' @export
read_oddi <- function(folder, meta_data_group, meta_data_xover = NULL, meta_data_sex = NULL, trim_bad_probe = T){
  tryCatch({
    files <- list.files(folder)
    samples <- sub(files, pattern = ".xlsx", replacement = "")
    if(is.null(meta_data_sex)){meta_data_sex <- tibble("Unknown" = samples)}
    if(is.null(meta_data_xover)){meta_data_xover <- tibble("NA" = samples)}
    temp <- list()
    for(i in files){
      temp[[i]] <- readxl::read_excel(paste0(folder, i))
      temp[[i]]$Mouse <- i
      colnames(temp[[i]])[2] <- "DegC"
      temp[[i]]$Group <- colnames(meta_data_group[grep(sub(i, pattern = ".xlsx", replacement = ""),meta_data_group)])
      temp[[i]]$Xover <- colnames(meta_data_xover[grep(sub(i, pattern = ".xlsx", replacement = ""),meta_data_xover)])
      temp[[i]]$Sex <- colnames(meta_data_sex[grep(sub(i, pattern = ".xlsx", replacement = ""),meta_data_sex)])
    }
    ret <- bind_rows(temp)
    ret$Mouse <- sub(ret$Mouse, pattern = ".xlsx", replacement = "") %>% factor()
    ret$Group <-ret$Group %>% factor()
    ret$Xover <-ret$Xover %>% factor()
    ret$Sex <-ret$Sex %>% factor()
    ret <- mutate(ret, Time = as.POSIXct("1900-01-01 00:00:00") + ((`Date & Time`)*(86400)))
    ret$Time <- lubridate::round_date(ret$Time, unit = "1 minute")
    ret <- select(ret, -c(`Date & Time`))
    if(trim_bad_probe == T){
      ret <- trim_bad_probe(ret)
    }
    return(ret)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#' Read tidy telemetry tibbles and export them as a multi page excel file.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param filename filename and path for saving excel file
#' @param return_list logical: if true, returns multi page object as a list
#' @return list or null
#' @export
export_telem <- function(tidy_telem, filename, return_list = F){
  tryCatch({
    tidy_telem$Time <- as.character(tidy_telem$Time)
    tlist <- list()
    for(i in levels(tidy_telem$Group)){
      tlist[[paste0("Deg.C, ", i)]] <- tidy_telem %>%
        filter(Group == i) %>%
        select(-Counts, -Xover, -Sex, -Group) %>%
        pivot_wider(names_from = Mouse, values_from = DegC)
      tlist[[paste0("Counts, ", i)]] <- tidy_telem %>%
        filter(Group == i) %>%
        select(-DegC, -Xover, -Sex, -Group) %>%
        pivot_wider(names_from = Mouse, values_from = Counts)
    }
    writexl::write_xlsx(tlist, path = filename)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  if(return_list ==T){return(tlist)}
}
#' Graph tidy telemetry data.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param telem_var telemetry measurement to plot.
#' @param one_day_avg collapse data into a 1 day average?
#' @param group_by variable to group by. Either "Group" or "Mouse"
#' @param tx_time vector of times (in posixct format) where to draw vertical dashed lines
#' @param se_ribbon logical: display standard error ribbon?
#' @return plot
#' @export
graph_telem <- function(tidy_telem, telem_var = "DegC", one_day_avg = F, group_by = "Mouse", tx_time = NULL, se_ribbon = T){
  tryCatch({
    if(is.POSIXct(tidy_telem$Time)){
      if(one_day_avg == F){
        tidy_telem$DST <- lubridate::dst(tidy_telem$Time)
        split1 <- filter(tidy_telem, DST == T)
        split2 <- filter(tidy_telem, DST == F)
        rect_right1 <- split1$Time
        lubridate::hour(rect_right1) <- 0
        lubridate::minute(rect_right1) <- 0
        rect_right1 <- rect_right1 %>% unique() + lubridate::hours(6)
        rect_right2 <- split2$Time
        lubridate::hour(rect_right2) <- 0
        lubridate::minute(rect_right2) <- 0
        rect_right2 <- rect_right2 %>% unique() + lubridate::hours(7)
        rect_right <- c(rect_right1, rect_right2)
        rect_right <- rect_right[rect_right %within% interval(first_time(tidy_telem$Time), last_time(tidy_telem$Time))]
        rectangles <- data.frame(
          xmin = rect_right - lubridate::hours(12),
          xmax = rect_right,
          ymin = -Inf,
          ymax = Inf)
        time1 <- c(first_time(grep("12:00:00", x = tidy_telem[["Time"]], value = T)), first_time(grep("00:00:00", x = tidy_telem[["Time"]], value = T))) %>%
          as.POSIXct() %>%
          first_time()
        time2 <- c(last_time(grep("12:00:00", x = tidy_telem[["Time"]], value = T)), last_time(grep("00:00:00", x = tidy_telem[["Time"]], value = T))) %>%
          as.POSIXct() %>%
          last_time()
        p1 <- ggplot(data = tidy_telem, aes(x = tidy_telem[["Time"]], y = tidy_telem[[telem_var]], color = tidy_telem[[group_by]], group = tidy_telem[[group_by]])) +
          theme_classic() +
          geom_rect(inherit.aes = F, data = rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    fill = 'gray80', alpha = 0.8) +
          stat_summary(geom = "line", fun.y = mean) +
          xlab("Time") +
          ylab(telem_var) +
          scale_color_discrete(name = group_by) +
          scale_x_datetime(breaks = seq(time1, time2, by = 43200)) +
          theme(axis.text.x = element_text(angle = -90)) +
          labs(color = "Group")

      } else {
        tdiff <- last_time(tidy_telem[["Time"]]) - first_time(tidy_telem[["Time"]])
        tdiff <- tdiff %>%
          round(2) %>%
          as.numeric()
        tidy_telem <- collapse_telem(tidy_telem)
        last.break <- tidy_telem$Time %>% last() %>% as.character()
        p1 <- ggplot(data = tidy_telem, aes(x = tidy_telem[["Time"]], y = tidy_telem[[telem_var]], color = tidy_telem[[group_by]], group = tidy_telem[[group_by]])) +
          theme_classic() +
          geom_rect(aes(xmin = "12:00:00",
                        xmax = last.break,
                        ymin = -Inf, ymax = Inf), color = "white", fill = "lightgrey") +
          stat_summary(geom = "line", fun.y = mean) +
          xlab("Zeitgeber Time") +
          ylab(telem_var) +
          scale_x_discrete(breaks = c("00:00:00", "06:00:00", "12:00:00", "18:00:00", last.break)) +
          scale_color_discrete(name = group_by) +
          theme(axis.text.x = element_text(angle = -90)) +
          ggtitle(paste0(tdiff, " day mean"))  +
          labs(color = "Group")
      }
      if(!is.null(tx_time)){
        p1 <- p1 +
          geom_vline(xintercept = as.POSIXct(tx_time), linetype = "dashed")
      }
      if(se_ribbon == T){
        p1 <- p1 +
          stat_summary(geom="ribbon", fun.data = mean_se, aes(fill = tidy_telem[[group_by]]), alpha = .5, color = NA, show.legend = F)
      }
    } else {
      last.break <- tidy_telem$Time %>% last() %>% as.character()
      p1 <- ggplot(data = tidy_telem, aes(x = tidy_telem[["Time"]], y = tidy_telem[[telem_var]], color = tidy_telem[[group_by]], group = tidy_telem[[group_by]])) +
        theme_classic() +
        geom_rect(aes(xmin = "12:00:00",
                      xmax = last.break,
                      ymin = -Inf, ymax = Inf), color = "white", fill = "lightgrey") +
        stat_summary(geom = "line", fun.y = mean) +
        xlab("Zeitgeber Time") +
        ylab(telem_var) +
        scale_x_discrete(breaks = c("00:00:00", "06:00:00", "12:00:00", "18:00:00", last.break)) +
        scale_color_discrete(name = group_by) +
        theme(axis.text.x = element_text(angle = -90)) +
        labs(color = "Group")
      if(!is.null(tx_time)){
        p1 <- p1 +
          geom_vline(xintercept = as.POSIXct(tx_time), linetype = "dashed")
      }
      if(se_ribbon == T){
        p1 <- p1 +
          stat_summary(geom="ribbon", fun.data = mean_se, aes(fill = tidy_telem[[group_by]]), alpha = .5, color = NA, show.legend = F)
      }
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(p1)
}
#' Trim tidy telemetry data.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param start_time start time for new output file
#' @param end_time end time for new output file
#' @param keep_groups optionally keep selected groups as defined in meta data. Default (NULL) is to keep all.
#' @param keep_mice optionally keep defined mice. Default (NULL) is to keep all.
#' @return a tidy telemetry tibble
#' @export
trim_telem <- function(tidy_telem, start_time = NULL, end_time = NULL, keep_groups = NULL, keep_mice = NULL, cut_mice = NULL){
  tryCatch({
    if(!is.null(start_time)){
      tidy_telem <- tidy_telem %>%
        filter(Time >= start_time)
    }
    if(!is.null(end_time)){
      tidy_telem <- tidy_telem %>%
        filter(Time <= end_time)
    }
    if(!is.null(keep_groups)){
      tidy_telem <- tidy_telem %>%
        filter(Group %in% keep_groups)
    }
    if(!is.null(keep_mice)){
      tidy_telem <- tidy_telem %>%
        filter(Mouse %in% keep_mice)
    }
    if(!is.null(cut_mice)){
      tidy_telem <- tidy_telem %>%
        filter(!Mouse %in% cut_mice)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(tidy_telem)
}
#' linear mixed effects model of telemetry data.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param telem_var telemetry measurement to analyze
#' @param collapse_first logical. collapse to a 1 day average before analysis? recommended.
#' @return list containing anova and models built
#' @export
lme_telem <- function(tidy_telem, telem_var = "DegC", collapse_first = T){
  tryCatch({
    if(collapse_first == T){
      tidy_telem <- tidy_telem %>% collapse_telem()
    }
    basemjev <<- formula(paste0(telem_var," ~ 1"))
    groupmjev <<- formula(paste0(telem_var," ~ Group"))
    timemjev <<- formula(paste0(telem_var," ~ Group + Time"))
    intermjev <<- formula(paste0(telem_var," ~ Group * Time"))
    baseline <- nlme::lme(basemjev, random = ~1 | Mouse/Group, data = tidy_telem, method = "ML", na.action = na.omit)
    Group <- nlme::lme(groupmjev, random = ~1 | Mouse/Group, data = tidy_telem, method = "ML", na.action = na.omit)
    Time <- nlme::lme(timemjev, random = ~1 | Mouse/Group, data = tidy_telem, method = "ML", na.action = na.omit)
    GroupXTime <- nlme::lme(intermjev, random = ~1 | Mouse/Group, data = tidy_telem, method = "ML", na.action = na.omit)
    statlist <- list()
    statlist[["AOV"]] <- anova(baseline, Group, Time, GroupXTime)
    statlist[["GroupSummary"]] <- summary(Group)
    statlist[["TimeSummary"]] <- summary(Time)
    statlist[["InteractionSummary"]] <- summary(GroupXTime)
    rm(basemjev, groupmjev, timemjev, intermjev, pos = ".GlobalEnv")
    return(statlist)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#' multiple comparisons of telemetry data.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param formula a fourmula (unquoted) for comparisons. E.g. DegC ~ Group or Counts ~ Group
#' @param collapse_first optionally collapse to one day before multiple comparisons
#' @param p_adjust_method method for adjusting p val. one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @return plot or table of multiple comparisons
#' @export
multcomp_telem <- function(tidy_telem, formula, plot_or_table = "plot", collapse_first = T, p_adjust_method = "BH"){
  tryCatch({
    if(collapse_first == T){
      tidy_telem <- tidy_telem %>% collapse_telem()
    }
    multcomp <- tidy_telem %>%
      group_by(Time) %>%
      do(broom::tidy(t.test(formula, data = .))) %>%
      ungroup %>%
      mutate(padj = p.adjust(p.value, method = p_adjust_method))
    multcomp$stat_summary <- sapply(X = multcomp$padj, FUN = MCsumm)
    if(plot_or_table == "plot"){
      p1 <- ggplot(data = multcomp, aes(x = Time, y = -log10(padj))) +
        geom_point() +
        theme_classic() +
        geom_hline(aes(yintercept = 1.30103, color = "p = .05"), linetype = "dashed") +
        geom_hline(aes(yintercept = 2, color = "p = .01"), linetype = "dashed") +
        geom_hline(aes(yintercept = 3, color = "p = .001"), linetype = "dashed") +
        geom_hline(aes(yintercept = 4, color = "p = .0001"), linetype = "dashed") +
        theme(axis.text.x = element_text(angle = -90)) +
        labs(color = "p values")

      if(!is.POSIXct(tidy_telem$Time)){
        p1 <- p1 + scale_x_discrete(breaks = c("06:00:00", "12:00:00", "18:00:00", "00:00:00", "23:55:00")) +
          xlab("Zeitgeber Time")
      } else {
        first_time <- function(time_vector){range(time_vector)[1]}
        last_time <- function(time_vector){range(time_vector)[2]}
        time1 <- c(first_time(grep("12:00:00", x = tidy_telem[["Time"]], value = T)), first_time(grep("00:00:00", x = tidy_telem[["Time"]], value = T))) %>%
          as.POSIXct() %>%
          first_time()
        time2 <- c(last_time(grep("12:00:00", x = tidy_telem[["Time"]], value = T)), last_time(grep("00:00:00", x = tidy_telem[["Time"]], value = T))) %>%
          as.POSIXct() %>%
          last_time()
        p1 <- p1 + scale_x_datetime(breaks = seq(time1, time2, by = 43200))
      }

    }
    else {
      p1 <- multcomp
    }
    return(p1)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#' return asterisks from p values for easy reading
#'
#' @param pval a p value
#' @return an asterisk code indicating significance
MCsumm <- function(pval){
  if(pval >=.05){x <- "NS"}
  if(pval < .05 & pval >= .01){x <- "*"}
  if(pval < .01 & pval >= .001){x <- "**"}
  if(pval < .001 & pval >= .0001){x <- "***"}
  if(pval < .0001){x <- "****"}
  return(x)
}
#' average telemetry data over user defined windows. very useful when data appear very noisy/oversampled
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param output_window what should final output increment be?
#' @return a time fuzzed tidy telemetry tibble
#' @export
fuzz_telem <-function(tidy_telem, output_window){
  tryCatch({
    tidy_telem$Time <- lubridate::round_date(tidy_telem$Time, unit = output_window)
    if("Counts" %in% colnames(tidy_telem)){
      tidy_telem <- tidy_telem %>%
        group_by_at(setdiff(names(tidy_telem), c("DegC", "Counts"))) %>%
        summarise(Counts = mean(Counts, na.rm = T), DegC = mean(DegC, na.rm = T))
    }else{
      tidy_telem <- tidy_telem %>%
        group_by_at(setdiff(names(tidy_telem), c("DegC"))) %>%
        summarise(DegC = mean(DegC, na.rm = T))
    }
    tidy_telem <- ungroup(tidy_telem)
    return(tidy_telem)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#' automatically detect and drop bad probes once they stop working. very useful for probes attached to mouse tails as they can fall off.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param DegC_threshold what temperature in C is considered a spuriously low reading?
#' @param over_how_long temp values should be below the threshold over how big of a time window before probe is considered bad?
#' @return a tidy telemetry tibble with bad probes removed starting a the point that they meet "bad" criteria
#' @export
trim_bad_probe <- function(tidy_telem, DegC_threshold = 23, over_how_long = "4 hours"){
  tryCatch({
    fuzzed <- fuzz_telem(tidy_telem, output_window = over_how_long)
    for(i in levels(tidy_telem$Mouse)){
      this.ms.fuzzed <- fuzzed %>%
        filter(Mouse == i)
      if(min(this.ms.fuzzed$DegC) < DegC_threshold) {
        first.time <- first(this.ms.fuzzed$Time)
        last.time.ms <- this.ms.fuzzed %>%
          filter(DegC < DegC_threshold)
        last.time <- last.time.ms$Time %>%
          first() %>%
          - lubridate::hours(sub(" hours", "", over_how_long) %>% as.numeric())
        this.ms.fixed <- tidy_telem %>%
          filter(Mouse == i) %>%
          trim_telem(start_time = first.time, end_time = last.time)
        tidy_telem <- tidy_telem %>%
          filter(Mouse != i) %>%
          bind_rows(this.ms.fixed)
      }}
    return(tidy_telem)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#' collapse all days in a tidy telem tibble, averaging read values
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @return a tidy telemetry tibble of length 1d. Note time is now of class factor
#' @export
collapse_telem <- function(tidy_telem){
  tryCatch({
    if(is.POSIXct(tidy_telem$Time)){
      tidy_telem <- zeitgeber_time(tidy_telem)
      t2 <- as.POSIXct("2020-01-11 23:55:00 PST")
      t3 <- as.POSIXct("2020-01-11 00:00:00 PST")
      levs <- c(seq(t3, t2, by=300)) %>%
        format("%H:%M:%S")
      tidy_telem$Time <- tidy_telem$Time %>%
        format("%H:%M:%S") %>%
        factor(levels = levs)
      if("Counts" %in% colnames(tidy_telem)){
        tidy_telem <- tidy_telem %>%
          group_by_at(setdiff(names(tidy_telem), c("DegC", "Counts"))) %>%
          summarise(Counts = mean(Counts, na.rm = T), DegC = mean(DegC, na.rm = T))
      }else{
        tidy_telem <- tidy_telem %>%
          group_by_at(setdiff(names(tidy_telem), c("DegC"))) %>%
          summarise(DegC = mean(DegC, na.rm = T))
      }
      tidy_telem <- ungroup(tidy_telem)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(tidy_telem)
}
#' convert time to zeitgeber (lightgiver) time: 0 = lights on, 12 = lights off.
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param undo logical. if true, converts ZT to clock time.
#' @param lights_on_DST at what time do the lights in the room come on, during DST. This algorithm assumes that the room light schedule does not adjust for DST
#' @return a tidy telemetry tibble converted to zeitgeber time.
#' @export
zeitgeber_time <- function(tidy_telem, undo = F, lights_on_DST = 8){
  tryCatch({
    if(undo == F){
      tidy_telem$DST <- lubridate::dst(tidy_telem$Time)
      split1 <- filter(tidy_telem, DST == T)
      split2 <- filter(tidy_telem, DST == F)
      split1$Time <- split1$Time %>% - lubridate::hours(lights_on_DST)
      split2$Time <- split2$Time %>% - lubridate::hours(lights_on_DST - 1)
      ret <- bind_rows(split1, split2) %>%
        select(-DST)
    } else {
      tidy_telem$DST <- lubridate::dst(tidy_telem$Time)
      split1 <- filter(tidy_telem, DST == T)
      split2 <- filter(tidy_telem, DST == F)
      split1$Time <- split1$Time %>% + lubridate::hours(lights_on_DST)
      split2$Time <- split2$Time %>% + lubridate::hours(lights_on_DST - 1)
      ret <- bind_rows(split1, split2) %>%
        select(-DST)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(ret)
}
#' add a logical column to tidy telem tibble indicating if lights are on (true)
#'
#' @param tidy_telem a tidy telemetry tibble, as produced by read_starr or read_oddi
#' @param lights_on_DST at what time do the lights in the room come on, during DST. This algorithm assumes that the room light schedule does not adjust for DST
#' @return a tidy telemetry tibble with logical column "DST" added. Time is still clock time.
light_dark <- function(tidy_telem, lights_on_DST = 8){
  tidy_telem <- tidy_telem %>%
    zeitgeber_time(lights_on_DST = lights_on_DST)
  tidy_telem$Light <- lubridate::am(tidy_telem$Time)
  tidy_telem <- tidy_telem %>%
    zeitgeber_time(undo = T, lights_on_DST = lights_on_DST)
  return(tidy_telem)
}
#' return the earliest time in a time vector
#'
#' @param time_vector a vector of times in class posixct or anything that works with range()
#' @return the earliest time in the vector
first_time <- function(time_vector){range(time_vector)[1]}
#' return the latest time in a time vector
#'
#' @param time_vector a vector of times in class posixct or anything that works with range()
#' @return the latest time in the vector
last_time <- function(time_vector){range(time_vector)[2]}
