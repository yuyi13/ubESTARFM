# Script: test-predict-kernel.R
# Objective: Verify R prediction kernels against independent published equations.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Small trained models and synthetic coarse target matrices.
# Outputs: testthat assertions for numerical prediction behavior.
# Usage: Rscript -e 'library(ubestarfm); testthat::test_file("tests/testthat/test-predict-kernel.R")'
# Dependencies: R packages testthat and ubestarfm.

test_that("detailed prediction matches the published equations", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 0.5
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
  actual <- ubestarfm:::ubestarfm_predict_arrays(
    model,
    list(target),
    c(-Inf, Inf)
  )[[1L]]

  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})

test_that("fallback prediction matches the published equations", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 0.5
  target[2:4, 2:4] <- NA_real_
  target[3, 3] <- fixture$coarse_1[3, 3] + 0.5
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
  actual <- ubestarfm:::ubestarfm_predict_arrays(
    model,
    list(target),
    c(-Inf, Inf)
  )[[1L]]

  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})

test_that("target missing values filter and renormalize candidates", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 0.5
  target[4, 4] <- NA_real_
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
  actual <- ubestarfm:::ubestarfm_predict_arrays(
    model,
    list(target),
    c(-Inf, Inf)
  )[[1L]]

  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})

test_that("single and batch prediction are identical", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  first <- fixture$coarse_1 + 0.5
  second <- fixture$coarse_1 + 1
  batch <- ubestarfm:::ubestarfm_predict_arrays(
    model,
    list(first, second),
    c(-Inf, Inf)
  )

  expect_equal(
    batch[[1L]],
    ubestarfm:::ubestarfm_predict_arrays(
      model,
      list(first),
      c(-Inf, Inf)
    )[[1L]]
  )
})

test_that("baseline and zero-bias modes match independent oracles", {
  fixture <- small_reference_fixture()
  target <- fixture$coarse_1 + 0.5
  for (method in c("baseline", "zero_bias")) {
    model <- do.call(
      ubestarfm:::ubestarfm_train_arrays,
      c(
        fixture,
        list(
          window_radius = 1L,
          patch_size = 5L,
          method = method
        )
      )
    )
    expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
    actual <- ubestarfm:::ubestarfm_predict_arrays(
      model,
      list(target),
      c(-Inf, Inf)
    )[[1L]]

    expect_equal(actual[3, 3], expected, tolerance = 1e-10)
  }
})

test_that("value-range correction uses corrected fine candidates", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 1000
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(0, 1))
  actual <- ubestarfm:::ubestarfm_predict_arrays(
    model,
    list(target),
    c(0, 1)
  )[[1L]]

  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})
