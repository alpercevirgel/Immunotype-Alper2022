
### %in& negate
`%!in%` = Negate(`%in%`)

### cor.mtest
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# colorblind friendly colors for clusters
cols_cluster <- c("1"= "#77AADD", "2"= "#99DDFF",
                  "3"= "#44BB99", "4"= "#BBCC33",
                  "5"= "#AAAA00", "6"= "#EEDD88",
                  "7"= "#EE8866", "8"= "#FFAABB", 
                  "9"= "#DDDDDD")

# colorblind friendly colors for age groups
cols_agegroup <- c(`25-49` = "#F0E442", `50-64` = "#85C0F9", `65-98` = "#F5793A")

#ggplot functions
boxplot_agegroups <- function()
{
  p <- ggplot(df_loop, aes(x=Xaxis, y=Yaxis)) +
    geom_boxplot(aes(fill=Groups), alpha=0.9,outlier.size=0,outlier.colour="white") +
    geom_jitter(aes(fill=Groups), alpha=0.4, width = 0.3, shape=21,size=1) +
    scale_fill_manual("Age groups", values=cols_agegroup) +
    theme_classic()+
    stat_summary(aes(y=Yaxis, x=Xaxis),size=0.2)
}

boxplot_agegroups_log <- function()
{
  p <- ggplot(df_loop, aes(x=Xaxis, y=Yaxis)) +
    geom_boxplot(aes(fill=Groups), alpha=0.9,outlier.size=0,outlier.colour="white") +
    geom_jitter(aes(fill=Groups), alpha=0.4, width = 0.3, shape=21,size=1) +
    scale_fill_manual("Age groups", values=cols_agegroup) +
    theme_classic()+
    stat_summary(aes(y=Yaxis, x=Xaxis),size=0.2)+
    scale_y_continuous(trans = 'log10') + annotation_logticks(sides="l")  
  
}

boxplot_cluster <- function()
{
  p <- ggplot(df_loop, aes(x=Xaxis, y=Yaxis)) +
    geom_boxplot(aes(fill=Groups), alpha=0.9,outlier.size=0,outlier.colour="white") +
    geom_jitter(aes(fill=Groups), alpha=0.4, width = 0.3, shape=21,size=1) +
    scale_fill_manual("Clusters", values=cols_cluster) +
    theme_classic()+
    stat_summary(aes(y=Yaxis, x=Xaxis),size=0.2)
}


boxplot_cluster_log <- function()
{
  p <- ggplot(df_loop, aes(x=Xaxis, y=Yaxis)) +
    geom_boxplot(aes(fill=Groups), alpha=0.9,outlier.size=0,outlier.colour="white") +
    geom_jitter(aes(fill=Groups), alpha=0.4, width = 0.3, shape=21,size=1) +
    scale_fill_manual("Clusters", values=cols_cluster) +
    theme_classic()+
    stat_summary(aes(y=Yaxis, x=Xaxis),size=0.2)+
    scale_y_continuous(trans = 'log10') + annotation_logticks(sides="l")  
}

### Plot dynamics function Day 0, Day 1, Day 2, Day 7

#' Plot Dynamics
#'
#' @name plotDynamics
#'
#' @details From a data.frame plots dynamics.
#'
#' @param input Data frame (wide format) with sample_identifier, subject_identifier,
#'          Timepoint, Age_group and trucount measurements. 
#'
#' @param variables a vector of variables to plot
#' @param save.dir path to save the out pdfs example "test_pdf/"
#' @param fig.heigth heigth of the figure 
#' @param fig.width width of the figure
#' @param color.vars Colors for group
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' library(reshape2)
#' # df_kinetics <- read_excel("df_kinetics.xlsx", 
#'                                        na = "n/a")
#' # variables.to.plot <- colnames(df_kinetics)[160:168]
#' 
#' # pick colors https://coolors.co/palettes/trending
#' # group.colors <- c( `25-49` = "#003049", `50-64` = "#118ab2", `65-98` = "#5f0f40")
#' # create a directory to store outputs
#' # dir.create("test_pdf")
#' plotDynamics(df_kinetics)
#' 
#' @return one or several pdfs with plots
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{Sudarshan
#' }}
#'
#' @export
NULL

group.colors_cluster <- c("1"= "#77AADD", "2"= "#99DDFF",
                          "3"= "#44BB99", "4"= "#BBCC33",
                          "5"= "#AAAA00", "6"= "#EEDD88",
                          "7"= "#EE8866", "8"= "#FFAABB", 
                          "9"= "#DDDDDD")

group.colors_agegroup <- c(`25_49` = "#F0E442", `50_64` = "#85C0F9", `65_98` = "#F5793A")


plotDynamics_cluster_ABC <- function(input=df_kinetics, 
                                     variables= variables.to.plot,
                                     save.dir = "results/figures/cluster_postvac_kinetics/",
                                     fig.heigth=4,
                                     fig.width = 6,
                                     color.vars = group.colors_cluster) {
  
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  
  first_df <- input %>% 
    dplyr::select(sample_identifier, 
                  #subject_identifier,
                  timepoint, cluster_number, 
                  !!!syms(variables)) %>% 
    reshape2::melt() %>% 
    mutate(Day= ifelse(timepoint=="A", 0, ifelse(timepoint =="B", 2, ifelse(timepoint =="C", 7, NA)))) 
  
  for(i in variables) {
    
    sub_df <- subset(first_df, variable == i)
    
    plasma_1 <- ggplot(sub_df, aes(x=as.numeric(Day), y=value, fill=cluster_number)) +
      stat_summary(geom="ribbon", alpha=0.6) +
      #annotate("label", x = 0, y = 0.5, label = "Day of vaccination")+
      geom_vline(xintercept = 0, lty=3) +
      geom_jitter(shape=21, size=2, alpha=0.2, width = 0.3,stroke = 0) +
      scale_fill_manual(values=color.vars) +
      #theme_biome_utils() +
      ylab(i) +
      xlab("Time (days)") +
      ggtitle(paste0("Dynamics of ", i)) +
      scale_y_log10()+
      theme_bw()
    
    filename_save <- paste0("Influenza_Dynamics_cluster_ABC_",i, ".pdf", sep="")
    save_path <- paste0(save.dir, filename_save, sep="")
    ggsave(save_path, height = fig.heigth, width = fig.width)
    
  }
}

plotDynamics_agegroup_ABC <- function(input=df_kinetics, 
                                      variables= variables.to.plot,
                                      save.dir = "results/figures/age_groups_postvac_kinetics/",
                                      fig.heigth=4,
                                      fig.width = 6,
                                      color.vars = group.colors_agegroup) {
  
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  
  first_df <- input %>% 
    dplyr::select(sample_identifier, 
                  #subject_identifier,
                  timepoint, age_group, 
                  !!!syms(variables)) %>% 
    reshape2::melt() %>% 
    mutate(Day= ifelse(timepoint=="A", 0, ifelse(timepoint =="B", 2, ifelse(timepoint =="C", 7,NA)))) 
  
  for(i in variables) {
    
    sub_df <- subset(first_df, variable == i)
    
    plasma_1 <- ggplot(sub_df, aes(x=as.numeric(Day), y=value, fill=age_group)) +
      stat_summary(geom="ribbon", alpha=0.6) +
      #annotate("label", x = 0, y = 0.5, label = "Day of vaccination")+
      geom_vline(xintercept = 0, lty=3) +
      geom_jitter(shape=21, size=2, alpha=0.2, width = 0.3,stroke = 0) +
      scale_fill_manual(values=color.vars) +
      #theme_biome_utils() +
      ylab(i) +
      xlab("Time (days)") +
      ggtitle(paste0("Dynamics of ", i)) +
      scale_y_log10()+
      theme_bw()
    
    filename_save <- paste0("Influenza_Dynamics_agegroup_ABC_",i, ".pdf", sep="")
    save_path <- paste0(save.dir, filename_save, sep="")
    ggsave(save_path, height = fig.heigth, width = fig.width)
    
  }
}





