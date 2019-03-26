rm(list= ls())
set.seed(19)

library(BAf)

N_SAMPLE <- 384
N_AB     <- 100
N_PLATE  <- N_SAMPLE / 96


##  BAf-like data generation ---------------------------------------------------
UPPER_LIMIT <- sample(19000:22000, 1)
LOWER_LIMIT <- sample(190:220, 1)

plate_effect  <- rnorm(N_PLATE,  sd= 1)
sample_effect <- rnorm(N_SAMPLE, sd= 1)
ab_effect     <- rnorm(N_AB, sd= 3)

#  EMPTY wells
i_empty <- c(10+c(0:(N_PLATE - 1))*97)
sample_effect[i_empty] <- rnorm(N_PLATE, mean= -9, sd= 0.1)

#  Replicated pools
i_pool <- sapply(c(0:(N_PLATE - 1)), function(ii) sample(96, 3) + ii * 96) %>%
    as.vector()
sample_effect[i_pool] <- sample_effect[i_pool] %>%
    mean(.) + rnorm(N_PLATE * 3, sd= 0.01)


#  negative control antibodies
ab_effect[N_AB:(N_AB - 1)] <- c(-8, -7)


X <- sapply(1:N_AB, function(ii) {           # per antibody
    # the effect of plate onto individual samples
    plate_effect_sample <- rep(plate_effect, each= N_SAMPLE / N_PLATE)
    ab_effect <- rep(ab_effect[ii], N_SAMPLE)
    ab_effect[i_empty] <- ab_effect[ii] * 0.1
    sample_effect + ab_effect + plate_effect_sample + rnorm(384, sd= 1)
}) %>%           # matrix generation
    { (. - min(.)) / (max(.) - min(.)) } %>%         #  0 <= X <= 1
    {. * log(UPPER_LIMIT / LOWER_LIMIT) + log(LOWER_LIMIT)} %>%   # 200 ~ 20000
    exp() %>%
    round()

colnames(X) <- paste0("T", 1:N_AB)
colnames(X)[N_AB:(N_AB - 2)] <- c("bare-bead", "rabbit IgG", "anti-human IgG")

#  positive control antibodies
X[, (N_AB - 2)] <- rnorm(N_SAMPLE, mean= UPPER_LIMIT + 10, sd= 500) %>%
    round()

plate_id <- seq_len(N_PLATE) %>%
    rep(each= 96) %>%
    factor()


##  Check distribution of values comparing plates ------------------------------

ask_do("check distribution", FUN= function() {
    geomean <- function(x) exp(mean(log(x)))
    tgm <- by(X, plate_id, apply, 2, geomean)
    
    opa <- par(mfrow= c(2, 3))
    apply(combn(N_PLATE, 2), 2, function(ea) {
        plot(tgm[[ea[1]]], tgm[[ea[2]]], 
             xlab= names(tgm)[ea[1]], ylab= names(tgm)[ea[2]], 
             log= "xy", asp= 1)
        abline(0, 1, col= "cadetblue")
    })
    par(opa)
    print(X[c(1:10, 97:106), 1:10])
    boxplot(t(X), log= "y")
})


##  Sample information ---------------------------------------------------------

change_half <- function(x, skip, value) {
    n <- length(x)
    x[sample(c(1:n)[-skip], (n - length(skip)) / 2)] <- value
    factor(x)
}

sInfo <- tibble("id"   = paste0("S", 1:N_SAMPLE), 
                "cohort"= rep("TEST", N_SAMPLE),
                "plate"= plate_id,
                "pos"  = rep(1:96, N_PLATE),
                "age"  = round(rbeta(N_SAMPLE, 2, 2) * 10 + 45),
                "sex"  = rep("female", N_SAMPLE) %>%
                    change_half(c(i_empty, i_pool), "male"),
                "dis"  = rep("cancer", N_SAMPLE) %>%
                    change_half(c(i_empty, i_pool), "no-cancer")
)
sInfo
summary(sInfo)

var_to_na <- c("age", "sex", "dis")

sInfo$id[i_empty] <- "EMPTY_0001"
sInfo$cohort[i_empty] <- "EMPTY"
sInfo[i_empty, var_to_na] <- NA

sInfo$id[i_pool]  <- "MIX_1_1004"
sInfo$cohort[i_pool] <- "MIX_1"
sInfo[i_pool,  var_to_na] <- NA


##  Binder info ------------------------------------------------------------------

binder <- tibble("assay"= rep(c("AY0", "AY1"), each= N_AB/2) %>% factor(),
                 "sba"  = rep(c("BA0", "BA1"), each= N_AB/2) %>% factor(),
                 "id" = colnames(X))



##  Assay-info for sample ------------------------------------------------------

assy_s <- list("fail_flag"= 
                 sapply(c("AY0", "AY1"), function(ii) {
                   rep("ok", N_SAMPLE) %>% {
                     .[sample(1:N_SAMPLE, 5)] <- "failed"
                     .
                   } %>%
                     factor()
                 }) %>% 
                 as.data.frame()
)

##  Assay-info for binder ------------------------------------------------------

assy_b <- matrix("ok", ncol= N_PLATE, nrow= N_AB) %>% {
  .[sample(1:N_AB, 2), ] <- "failed"
  .
} %>% 
  as.data.frame() %>% 
  `colnames<-`(levels(plate_id)) %>% 
  list("fail_flag"= .)

sba <- new("BAf", X, sinfo= sInfo, binder= binder, 
           sinfo_batch_i = "plate", binder_batch_i= "assay", 
           assay_sinfo = assy_s , assay_binder= assy_b)

devtools::use_data(sba, pkg= ".", overwrite= TRUE)

