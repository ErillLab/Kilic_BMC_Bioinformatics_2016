library(ggplot2)

df <- read.csv("../results_meme_on_pssm_searched_seeded.csv")
df <- df[tf=="LexA",]


p <- ggplot(df, aes(x=reference_target_dist,
                    y=target_inferred_dist,
                    )) +
            geom_point(alpha=0.5) +
            facet_wrap(~ transfer_method) +
            theme(legend.position="top")
p

ggsave(p, file="../euclidean_dist_old.png")

ggplot(df, aes(x=protein_distance,
               y=reference_target_dist)) +
    geom_point(alpha=0.5) +
    stat_smooth()
