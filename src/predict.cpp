#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

inline bool finite_at(const NumericMatrix& values, int row, int col) {
  return R_FINITE(values(row, col));
}

inline std::size_t row_major_index(int row, int col, int ncol) {
  return static_cast<std::size_t>(row) * ncol + col;
}

inline bool test_bit(
    const RawVector& masks,
    std::size_t base,
    int bit) {
  const std::size_t byte = base + static_cast<std::size_t>(bit / 8);
  const unsigned char value = static_cast<unsigned char>(masks[byte]);
  return (value & (1u << (bit % 8))) != 0;
}

inline double target_at(
    const NumericVector& targets,
    int target,
    int row,
    int col,
    int nrow,
    int ncol) {
  const std::size_t raster_size =
      static_cast<std::size_t>(nrow) * ncol;
  return targets[
      static_cast<std::size_t>(target) * raster_size +
      row_major_index(row, col, ncol)];
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix ubestarfm_predict_patch_cpp(
    const Rcpp::NumericMatrix& fine_1,
    const Rcpp::NumericMatrix& fine_2,
    const Rcpp::NumericMatrix& coarse_1,
    const Rcpp::NumericMatrix& coarse_2,
    const Rcpp::RawVector& candidate_masks,
    int mask_bytes,
    const Rcpp::NumericVector& targets,
    int target_count,
    int row_start,
    int row_end,
    int col_start,
    int col_end,
    int window_radius,
    std::string method,
    double value_min,
    double value_max) {
  const int nrow = fine_1.nrow();
  const int ncol = fine_1.ncol();
  const int side = 2 * window_radius + 1;
  const int patch_rows = row_end - row_start + 1;
  const int patch_cols = col_end - col_start + 1;
  const int patch_cells = patch_rows * patch_cols;
  NumericMatrix output(patch_cells, target_count);
  std::fill(output.begin(), output.end(), NA_REAL);

  int local_cell = 0;
  for (int row = row_start; row <= row_end; ++row) {
    for (int col = col_start; col <= col_end; ++col) {
      const bool center_reference_valid =
          finite_at(fine_1, row, col) &&
          finite_at(fine_2, row, col) &&
          finite_at(coarse_1, row, col) &&
          finite_at(coarse_2, row, col);
      if (!center_reference_valid) {
        ++local_cell;
        continue;
      }

      const std::size_t mask_base =
          row_major_index(row, col, ncol) * mask_bytes;
      const int window_row_start = std::max(0, row - window_radius);
      const int window_row_end = std::min(nrow - 1, row + window_radius);
      const int window_col_start = std::max(0, col - window_radius);
      const int window_col_end = std::min(ncol - 1, col + window_radius);

      for (int target = 0; target < target_count; ++target) {
        const double center_target = target_at(
            targets,
            target,
            row,
            col,
            nrow,
            ncol);
        if (!R_FINITE(center_target)) {
          continue;
        }

        double target_window_sum = 0.0;
        double coarse_1_window_sum = 0.0;
        double coarse_2_window_sum = 0.0;
        double fine_1_window_sum = 0.0;
        double fine_2_window_sum = 0.0;
        int window_count = 0;

        for (
            int window_col = window_col_start;
            window_col <= window_col_end;
            ++window_col) {
          for (
              int window_row = window_row_start;
              window_row <= window_row_end;
              ++window_row) {
            if (
                !finite_at(fine_1, window_row, window_col) ||
                !finite_at(fine_2, window_row, window_col) ||
                !finite_at(coarse_1, window_row, window_col) ||
                !finite_at(coarse_2, window_row, window_col)) {
              continue;
            }
            const double target_value = target_at(
                targets,
                target,
                window_row,
                window_col,
                nrow,
                ncol);
            if (!R_FINITE(target_value)) {
              continue;
            }

            target_window_sum += target_value;
            coarse_1_window_sum += coarse_1(window_row, window_col);
            coarse_2_window_sum += coarse_2(window_row, window_col);
            fine_1_window_sum += fine_1(window_row, window_col);
            fine_2_window_sum += fine_2(window_row, window_col);
            ++window_count;
          }
        }

        std::vector<int> candidate_rows;
        std::vector<int> candidate_cols;
        candidate_rows.reserve(side * side);
        candidate_cols.reserve(side * side);

        for (int delta_col = -window_radius; delta_col <= window_radius; ++delta_col) {
          for (int delta_row = -window_radius; delta_row <= window_radius; ++delta_row) {
            const int bit =
                (delta_col + window_radius) * side +
                (delta_row + window_radius);
            if (!test_bit(candidate_masks, mask_base, bit)) {
              continue;
            }

            const int candidate_row = row + delta_row;
            const int candidate_col = col + delta_col;
            if (
                candidate_row < 0 ||
                candidate_row >= nrow ||
                candidate_col < 0 ||
                candidate_col >= ncol) {
              continue;
            }
            const double candidate_target = target_at(
                targets,
                target,
                candidate_row,
                candidate_col,
                nrow,
                ncol);
            if (!R_FINITE(candidate_target)) {
              continue;
            }

            candidate_rows.push_back(candidate_row);
            candidate_cols.push_back(candidate_col);
          }
        }

        const int candidate_count =
            static_cast<int>(candidate_rows.size());
        if (candidate_count > 5) {
          std::vector<double> fine_candidates_1(candidate_count);
          std::vector<double> fine_candidates_2(candidate_count);
          std::vector<double> coarse_candidates_1(candidate_count);
          std::vector<double> coarse_candidates_2(candidate_count);
          std::vector<double> target_candidates(candidate_count);

          double fine_mean_1 = 0.0;
          double fine_mean_2 = 0.0;
          double coarse_mean_1 = 0.0;
          double coarse_mean_2 = 0.0;
          for (int candidate = 0; candidate < candidate_count; ++candidate) {
            const int candidate_row = candidate_rows[candidate];
            const int candidate_col = candidate_cols[candidate];
            fine_candidates_1[candidate] =
                fine_1(candidate_row, candidate_col);
            fine_candidates_2[candidate] =
                fine_2(candidate_row, candidate_col);
            coarse_candidates_1[candidate] =
                coarse_1(candidate_row, candidate_col);
            coarse_candidates_2[candidate] =
                coarse_2(candidate_row, candidate_col);
            target_candidates[candidate] = target_at(
                targets,
                target,
                candidate_row,
                candidate_col,
                nrow,
                ncol);
            fine_mean_1 += fine_candidates_1[candidate];
            fine_mean_2 += fine_candidates_2[candidate];
            coarse_mean_1 += coarse_candidates_1[candidate];
            coarse_mean_2 += coarse_candidates_2[candidate];
          }
          fine_mean_1 /= candidate_count;
          fine_mean_2 /= candidate_count;
          coarse_mean_1 /= candidate_count;
          coarse_mean_2 /= candidate_count;

          double fine_pixel_1 = fine_1(row, col);
          double fine_pixel_2 = fine_2(row, col);
          if (method == "zero_bias") {
            const double bias_1 = -fine_mean_1 + coarse_mean_1;
            const double bias_2 = -fine_mean_2 + coarse_mean_2;
            fine_pixel_1 += bias_1;
            fine_pixel_2 += bias_2;
            for (int candidate = 0; candidate < candidate_count; ++candidate) {
              fine_candidates_1[candidate] += bias_1;
              fine_candidates_2[candidate] += bias_2;
            }
          }

          std::vector<double> inverse_distance(candidate_count);
          double inverse_distance_sum = 0.0;
          for (int candidate = 0; candidate < candidate_count; ++candidate) {
            double spectral_distance =
                1.0 - 0.5 * (
                    std::abs(
                        (fine_candidates_1[candidate] -
                         coarse_candidates_1[candidate]) /
                        (fine_candidates_1[candidate] +
                         coarse_candidates_1[candidate])) +
                    std::abs(
                        (fine_candidates_2[candidate] -
                         coarse_candidates_2[candidate]) /
                        (fine_candidates_2[candidate] +
                         coarse_candidates_2[candidate])));
            if (spectral_distance > 1.0 || spectral_distance < -1.0) {
              spectral_distance = 0.5;
            }
            const double spatial_distance =
                1.0 + std::sqrt(
                    std::pow(col - candidate_cols[candidate], 2.0) +
                    std::pow(row - candidate_rows[candidate], 2.0)) /
                    window_radius;
            const double combined_distance =
                (1.0 - spectral_distance) * spatial_distance + 1e-7;
            inverse_distance[candidate] = 1.0 / combined_distance;
            inverse_distance_sum += inverse_distance[candidate];
          }

          const double mean_target =
              target_window_sum / window_count;
          const double mean_coarse_1 =
              coarse_1_window_sum / window_count;
          const double mean_coarse_2 =
              coarse_2_window_sum / window_count;
          const double temporal_difference_1 =
              std::abs(mean_target - mean_coarse_1) + 1e-10;
          const double temporal_difference_2 =
              std::abs(mean_target - mean_coarse_2) + 1e-10;
          const double inverse_temporal_1 =
              1.0 / temporal_difference_1;
          const double inverse_temporal_2 =
              1.0 / temporal_difference_2;
          const double temporal_weight_1 =
              inverse_temporal_1 /
              (inverse_temporal_1 + inverse_temporal_2);
          const double temporal_weight_2 = 1.0 - temporal_weight_1;

          double change_1 = 0.0;
          double change_2 = 0.0;
          double fallback_1 = 0.0;
          double fallback_2 = 0.0;
          for (int candidate = 0; candidate < candidate_count; ++candidate) {
            const double weight =
                inverse_distance[candidate] / inverse_distance_sum;
            change_1 += weight * (
                target_candidates[candidate] -
                coarse_candidates_1[candidate]);
            change_2 += weight * (
                target_candidates[candidate] -
                coarse_candidates_2[candidate]);
            fallback_1 += weight * fine_candidates_1[candidate];
            fallback_2 += weight * fine_candidates_2[candidate];
          }

          double prediction =
              temporal_weight_1 * (fine_pixel_1 + change_1) +
              temporal_weight_2 * (fine_pixel_2 + change_2);
          if (prediction <= value_min || prediction >= value_max) {
            prediction =
                temporal_weight_1 * fallback_1 +
                temporal_weight_2 * fallback_2;
          }
          output(local_cell, target) = prediction;
        } else {
          const double mean_target =
              target_window_sum / window_count;
          const double mean_coarse_1 =
              coarse_1_window_sum / window_count;
          const double mean_coarse_2 =
              coarse_2_window_sum / window_count;
          double fine_pixel_1 = fine_1(row, col);
          double fine_pixel_2 = fine_2(row, col);
          if (method == "zero_bias") {
            fine_pixel_1 =
                fine_pixel_1 -
                fine_1_window_sum / window_count +
                mean_coarse_1;
            fine_pixel_2 =
                fine_pixel_2 -
                fine_2_window_sum / window_count +
                mean_coarse_2;
          }

          const double temporal_difference_1 =
              mean_target - mean_coarse_1 + 1e-10;
          const double temporal_difference_2 =
              mean_target - mean_coarse_2 + 1e-10;
          const double inverse_temporal_1 =
              1.0 / std::abs(temporal_difference_1);
          const double inverse_temporal_2 =
              1.0 / std::abs(temporal_difference_2);
          const double temporal_weight_1 =
              inverse_temporal_1 /
              (inverse_temporal_1 + inverse_temporal_2);
          const double temporal_weight_2 = 1.0 - temporal_weight_1;
          output(local_cell, target) =
              temporal_weight_1 *
                  (fine_pixel_1 + temporal_difference_1) +
              temporal_weight_2 *
                  (fine_pixel_2 + temporal_difference_2);
        }
      }

      ++local_cell;
    }
  }

  return output;
}
