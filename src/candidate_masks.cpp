#include <Rcpp.h>
#include <cmath>
#include <cstddef>

using namespace Rcpp;

namespace {

inline std::size_t cell_index(int row, int col, int ncol) {
  return static_cast<std::size_t>(row) * ncol + col;
}

inline int local_bit_index(int delta_row, int delta_col, int radius) {
  const int side = 2 * radius + 1;
  return (delta_col + radius) * side + (delta_row + radius);
}

inline void set_bit(RawVector& masks, std::size_t base, int bit) {
  const std::size_t byte = base + static_cast<std::size_t>(bit / 8);
  const unsigned char value = static_cast<unsigned char>(masks[byte]);
  masks[byte] = static_cast<Rbyte>(value | (1u << (bit % 8)));
}

inline bool finite_at(const NumericMatrix& values, int row, int col) {
  return R_FINITE(values(row, col));
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List ubestarfm_train_patch_cpp(
    const Rcpp::NumericMatrix& fine_1,
    const Rcpp::NumericMatrix& fine_2,
    const Rcpp::NumericMatrix& coarse_1,
    const Rcpp::NumericMatrix& coarse_2,
    int row_start,
    int row_end,
    int col_start,
    int col_end,
    int window_radius,
    double threshold_1,
    double threshold_2) {
  const int nrow = fine_1.nrow();
  const int ncol = fine_1.ncol();
  const int side = 2 * window_radius + 1;
  const int mask_bytes = (side * side + 7) / 8;
  const int patch_rows = row_end - row_start + 1;
  const int patch_cols = col_end - col_start + 1;
  const int patch_cells = patch_rows * patch_cols;

  IntegerVector cell_ids(patch_cells);
  RawVector masks(static_cast<R_xlen_t>(patch_cells) * mask_bytes);

  int local_cell = 0;
  for (int row = row_start; row <= row_end; ++row) {
    for (int col = col_start; col <= col_end; ++col) {
      cell_ids[local_cell] = static_cast<int>(cell_index(row, col, ncol));
      const std::size_t base =
          static_cast<std::size_t>(local_cell) * mask_bytes;

      const bool center_valid =
          finite_at(fine_1, row, col) &&
          finite_at(fine_2, row, col) &&
          finite_at(coarse_1, row, col) &&
          finite_at(coarse_2, row, col);

      if (center_valid) {
        const int candidate_row_start = std::max(0, row - window_radius);
        const int candidate_row_end =
            std::min(nrow - 1, row + window_radius);
        const int candidate_col_start = std::max(0, col - window_radius);
        const int candidate_col_end =
            std::min(ncol - 1, col + window_radius);

        for (
            int candidate_col = candidate_col_start;
            candidate_col <= candidate_col_end;
            ++candidate_col) {
          for (
              int candidate_row = candidate_row_start;
              candidate_row <= candidate_row_end;
              ++candidate_row) {
            if (
                !finite_at(fine_1, candidate_row, candidate_col) ||
                !finite_at(fine_2, candidate_row, candidate_col) ||
                !finite_at(coarse_1, candidate_row, candidate_col) ||
                !finite_at(coarse_2, candidate_row, candidate_col)) {
              continue;
            }

            if (
                std::abs(
                    fine_1(candidate_row, candidate_col) -
                    fine_1(row, col)) < threshold_1 &&
                std::abs(
                    fine_2(candidate_row, candidate_col) -
                    fine_2(row, col)) < threshold_2) {
              const int bit = local_bit_index(
                  candidate_row - row,
                  candidate_col - col,
                  window_radius);
              set_bit(masks, base, bit);
            }
          }
        }
      }

      ++local_cell;
    }
  }

  return List::create(
      Named("cell_ids") = cell_ids,
      Named("masks") = masks);
}

// [[Rcpp::export]]
Rcpp::RawVector ubestarfm_assemble_masks_cpp(
    const Rcpp::List& patch_results,
    int total_cells,
    int mask_bytes) {
  RawVector assembled(static_cast<R_xlen_t>(total_cells) * mask_bytes);

  for (R_xlen_t patch_index = 0; patch_index < patch_results.size(); ++patch_index) {
    const List patch = patch_results[patch_index];
    const IntegerVector cell_ids = patch["cell_ids"];
    const RawVector masks = patch["masks"];

    for (R_xlen_t local_cell = 0; local_cell < cell_ids.size(); ++local_cell) {
      const std::size_t source_base =
          static_cast<std::size_t>(local_cell) * mask_bytes;
      const std::size_t destination_base =
          static_cast<std::size_t>(cell_ids[local_cell]) * mask_bytes;
      for (int byte = 0; byte < mask_bytes; ++byte) {
        assembled[destination_base + byte] = masks[source_base + byte];
      }
    }
  }

  return assembled;
}
