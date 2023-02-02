rank_fh <- function(serial) {
  position <- seq(length(serial))
  #ranking <- is.na(serial)
  df <- data.frame(serial, position)
  df_NA <- df %>% filter(is.na(serial))
  if (nrow(df_NA) > 0) {
    df_NA$ranking = NA
  }
  df_Inf <- df %>% filter(is.infinite(serial))
  if (nrow(df_Inf) > 0) {
    df_Inf$ranking = NA
  }
  df_pos <- df %>% filter(serial > 0)
  if (nrow(df_pos) > 0) {
    df_pos <- df_pos %>% arrange(serial)
    df_pos$ranking <- rank(-df_pos$serial, ties.method = 'min')
  }
  df_zero <- df %>% filter(serial==0)
  if (nrow(df_zero) > 0) {
    df_zero$ranking = 0
  }
  df_neg <- df %>% filter(serial < 0)
  if (nrow(df_neg) > 0) {
    df_neg <- df_neg %>% arrange(serial)
    df_neg$ranking <- rank(df_neg$serial, ties.method = 'min')
  }
  df <- rbind(df_NA,df_pos,df_zero,df_neg)
  df <- df[order(df$position),]
  return(df$ranking)
}