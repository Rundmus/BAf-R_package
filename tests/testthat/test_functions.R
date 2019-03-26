expect_dim <- function(obj, d) {
  stopifnot(is.numeric(d))
  eval.parent(substitute(expect_equal(dim(obj), d)))
}

data(sba)

context("Read functions")

test_that("read_LIMS_SBA_files function", {
  lims_B1 <- read_LIMS_binder_file("../SBA_files/SL0005_binder.tsv", 1)
  expect_s3_class(lims_B1$multiple_target, "factor")
  expect_equal(levels(lims_B1$multiple_target), c("0", "1"))
  expect_s3_class(lims_B1, "data.frame")
  
  expect_error(read_LIMS_SBA_files(path = "../SBA_files"))
  expect_error(read_LIMS_SBA_files(SLid = "SL0005"))
  
  lims <- read_LIMS_SBA_files(SLid = "SL0005", path = "../SBA_files", 
                              QC_plot= F)
  expect_s4_class(lims, "BAf")
  expect_dim(lims, c(288, 800))
  expect_equal(levels(factor(batch(lims, "b"))), c("AY098", "AY249"))
})

test_that("parse_FlexMAP3D_csv_2_list function", {
  rFM <- parse_FlexMAP3D_csv_2_list(
    "../SBA_files/130531_SBA65_LYM_FH_clyde_20130531_145913.csv.gz"
  )
  expect_type(rFM, "list")
  expect_named(rFM, 
               c('head', 'pos', 'median', 'net_mfi', 'count', 'avg_net_mfi', 
                 'units', 'per_bead_count', 'dilution_factor', 
                 'analysis_types', 'audit_logs', 'warnings/errors')
  )
})




context("Other functions")



test_that("apply_per_group functions", {
  expect_dim(apply_per_group(sba, median, na.rm= TRUE, by_s= c("plate", "sex")),
             c(8, 1))
})


test_that("clean up functions", {
  expect_output(failed_to_NA(sba, wise= "sinfo", by= "plate"))
  expect_output(failed_to_NA(sba, by= "plate"))
  expect_output(failed_to_NA(sba, wise= "b", by= "assay"))
  expect_output(failed_to_NA(sba))
  expect_output(replace_0(sba))
  
})

test_that("dividyBy & meanBy", {
  expect_is(divideBy(sba, "plate"), "list")
  expect_length(divideBy(sba, "plate"), 5)
  expect_is(meanBy(sba, "plate"), "matrix")
  expect_is(meanBy(sba, "dis", geometric = TRUE), "matrix")
})

test_that("Repeat related functions", {
  sba2 <- cbind(sba[, 1:5], sba[, 1:5] + 1)
  expect_s4_class(combine_repeats(sba2, rep(bid(sba)[1:5], 2)), "BAf")
})

test_that("Normalization functions", {
  expect_is(pqn(sba@.Data), "matrix")
  expect_is(pqn(sba@.Data)[, 1], "numeric")
  
  expect_s4_class(pqn(sba), "BAf")
})






context("Plot functions")

test_that("Fundamental plot_QC functions", {
  expect_is(plot_QC_sample_signal_boxplot(sba), "list")
  expect_is(plot_QC_binder_signal_boxplot(sba), "list")
  expect_is(plot_QC_repl_var(sba, sI(sba)$cohort == "MIX_1"), "list")
  expect_is(
    plot_QC_bead_count_boxplot("../SBA_files/SL0005_data_bead_count.tsv"), 
    "list"
  )
})

test_that("Plot_PC functions", {
  expect_is(plot_PC(sba), "prcomp")
  expect_is(
    plot_PC(sba, col_by= "plate", pch_by= "sex", legend_arg= list()),
    "prcomp"
  )
})
