expect_dim <- function(obj, d) {
  stopifnot(is.numeric(d))
  eval.parent(substitute(expect_equal(dim(obj), d)))
}

data(sba)
N_SAMPLE <- 384
N_BINDER <- 100

context("BAf-class")

test_that("Create a new BAf object", {
  expect_s4_class(BAf(), "BAf")
  expect_error(BAf(sinfo= tibble(id= paste0("s", 1:10))))
})

test_that("Class validity", {
  expect_true(validObject(sba))
  expect_is(sba, "BAf")
  expect_is(sba, "matrix")
  tmp <- sba
  tmp@sinfo$id[3] <- "err"
  expect_error(validObject(tmp))
  tmp <- sba
  tmp@binder$id <- "err"
  expect_error(validObject(tmp))
  expect_error(sba@batch_c <- "err")
})

test_that("key IDs", {
  expect_equal(s_b_switch("s"), "sinfo")
  expect_equal(s_b_switch("binder"), "binder")
  expect_error(s_b_switch("sa"))
  
  expect_is(sid(sba), "character")
  expect_is(sid(sba, exact= T), "character")
  expect_length(sid(sba), N_SAMPLE)
  expect_is(bid(sba), "character")
  expect_length(bid(sba), N_BINDER)
  tmp <- sba
  sid(tmp) <- rep("SID", N_SAMPLE)
  expect_equal(sid(tmp, exact= T)[1:3], c("SID", "SID*1", "SID*2"))
  tmp <- sba
  bid(tmp) <- rep("BID", N_BINDER)
  expect_equal(bid(tmp, exact= T)[1:3], c("BID", "BID*1", "BID*2"))
  # expect_equal(.key_id_to_index(sba, 10), 10)
  # expect_equal(.key_id_to_index(sba, 1:10), 1:10)
  # expect_equal(.key_id_to_index(sba, 10, wise= "binder"), 10)
  # expect_equal(.key_id_to_index(sba, c("T24", "T21"), wise= "binder"), c(24, 21))
  # expect_warning(.key_id_to_index(sba, "S24*1", exact= F))
})

test_that("batch ID acess", {
  
  expect_equal(batch(sba, "binder"), rep(c("AY0", "AY1"), each= 50))
  batch(sba, "binder") <- rep(c("AY10", "AY11"), each= ncol(sba)/2)
  expect_equal(batch(sba, "binder"), rep(c("AY10", "AY11"), each= 50))
  batch(sba, "binder")[1:50] <- "AY12"
  expect_equal(batch(sba, "binder"), rep(c("AY12", "AY11"), each= 50))
  expect_equal(batch(sba, "b"), rep(c("AY12", "AY11"), each= 50))
  expect_equal(batch(sba, "sinfo"), as.character(rep(1:4, each= 96)))
  batch(sba, "sinfo") <- rep(LETTERS[1:4], each= 96)
  expect_equal(batch(sba, "sinfo"), rep(LETTERS[1:4], each= 96))
  batch(sba, "sinfo")[1:96] <- "F"
  expect_equal(batch(sba, "sinfo"), rep(c("F", "B", "C", "D"), each= 96))
  expect_equal(batch(sba, "sample"), rep(c("F", "B", "C", "D"), each= 96))
  
  expect_equal(batch_colname(sba), "plate")
  expect_equal(batch_colname(sba, "binder"), "assay")
})

test_that("Arithmetic computation", {
  expect_s4_class(sba + 0.5, "BAf")
  expect_identical((sba + 0.5)@.Data, (sba@.Data + 0.5))
})


context("Basic functions")

test_that("data access", {
  expect_is(sX(sba), "matrix")
  expect_is(sI(sba), "data.frame")
  expect_is(sB(sba), "data.frame")
  expect_is(sA(sba)[[1]], "data.frame")
  expect_equal(sA(sba)[[1]], sA(sba, "sinfo")[[1]])
  expect_is(codebook(sba)$sinfo, "data.frame")
  expect_is(note(sba), "character")
})

test_that("Subset (or extract)", {
  expect_dim(sba[0, ], c(0, N_BINDER))
  expect_dim(sba[0, ][, 1:10], c(0, 10))
  expect_dim(sba[1:4, 3], c(4, 1))
  expect_equal(dim(sba[, c("T2", "T3")])[2L], 2)
})

test_that("cbind / rbind", {
  expect_dim(cbind(sba[, 1:3], sba[, 10:20]), c(384, 14))
  expect_dim(rbind(sba[1:5, ], sba[10:20, ]), c(16, 100))
})

test_that("aux_arg_handler.R", {
  expect_equal(length(index_grouped_by_cat(sba@sinfo, "plate")), 4)
  expect_equivalent(sapply(index_grouped_by_cat(sba@sinfo, "plate"), length), 
                    rep(96, 4))
  expect_identical(
    anyDuplicated(unlist(index_grouped_by_cat(sba@sinfo, "plate"))), 0L
  )
})


