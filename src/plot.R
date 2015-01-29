library(ggplot2)

df <- read.csv("../results.csv",
               colClasses=c("reference_target_dist"="numeric",
                   "target_inferred_dist"="numeric",
                   "protein_distance"="numeric"))
df <- df[df$transfer_method != "random",]

ggplot(df, aes(x=reference_target_dist,
               y=target_inferred_dist,
               color=transfer_method)) +
    geom_point(alpha=0.5) +
    stat_smooth() +
    facet_wrap(~ transfer_method)



ggplot(df, aes(x=protein_distance,
               y=reference_target_dist)) +
    geom_point(alpha=0.5) +
    stat_smooth()
