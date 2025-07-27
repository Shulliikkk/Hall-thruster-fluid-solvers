#ifndef Interp
#define Interp

#include <vector>
#include <cmath>

namespace interp {
  // Class for working with SLAE with a tridiagonal matrix
  template<typename T>
  class Tridiagonal_System {
  private:
      std::vector<T> a, c;
      std::vector<T> b, column;

  public:
      Tridiagonal_System(const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c, const std::vector<T>& column) : a(a), b(b), c(c), column(column) {}

      // Tridiagonal matrix algorithm
      std::vector<T> solve_tridig_mat_alg() {
          std::size_t N = b.size();
          std::vector<T> p(N - 1), q(N - 1);
          std::vector<T> x(N);

          p[0] = -c[0] / b[0];
          q[0] = column[0] / b[0];

          for (std::size_t i = 1; i < N - 1; i++) {
              p[i] = -c[i] / (a[i] * p[i - 1] + b[i]);
              q[i] = (column[i] - a[i] * q[i - 1]) / (a[i] * p[i - 1] + b[i]);
          }

          x[N - 1] = (column[N - 1] - a[N - 2] * q[N - 2]) / (p[N - 2] * a[N - 2] + b[N - 1]);

          for (std::size_t i = N - 2; i < -1; i--) {
              x[i] = p[i] * x[i + 1] + q[i];
          }

          return x;
      }
  };

  // Class of the natural cubic spline
  template<typename T>
  class Cubic_Spline {
  private:
      std::vector<T> A, B, C, D;
      std::vector<T> h, u, x;

      Tridiagonal_System<T> create_system(std::vector<T>& h, std::vector<T>& u) {
          std::size_t N = u.size();
          std::vector<T> a(N - 2), c(N - 2);
          std::vector<T> b(N - 1), column;

          T dvd_dif_1_0, dvd_dif_1_1;
          for (std::size_t i = 2; i < N; i++) {
              dvd_dif_1_0 = (u[i - 1] - u[i - 2]) / h[i - 2];
              dvd_dif_1_1 = (u[i] - u[i - 1]) / h[i - 1];
              column.push_back(6 * (dvd_dif_1_1 - dvd_dif_1_0) / (h[i - 1] + h[i - 2]));

              b[i - 2] = 2;

              if (i < N) {
                  c[i - 2] = h[i - 1] / (h[i - 2] + h[i - 1]);
              }
              if (i > 2) {
                  a[i - 2] = h[i - 2] / (h[i - 2] + h[i - 1]);
              }
          }

          return Tridiagonal_System<T>(a, b, c, column);
      }

  public:
      Cubic_Spline(const std::vector<T>& x, const std::vector<T>& y) : x(x), u(y), A(x.size()), B(x.size()), C(x.size()), D(x.size()), h(x.size()) {
          std::size_t N = x.size();

          for (std::size_t i = 1; i < N; i++) {
              h[i - 1] = x[i] - x[i - 1];
              A[i - 1] = u[i];
          }

          auto system = create_system(h, u);
          std::vector<T> C_0 = system.solve_tridig_mat_alg();

          for (std::size_t i = 0; i < N - 1; i++) {
              C[i] = C_0[i];
          }

          C[N - 1] = 0;
          D[0] = C[0] / h[0];

          for (std::size_t i = 1; i < N; i++) {
              D[i] = (C[i] - C[i - 1]) / h[i];
          }

          B[0] = A[0] / h[0] + C[0] / 2 * h[0] - D[0] / 6 * h[0] * h[0] - u[0] / h[0];

          for (std::size_t i = 1; i < N; i++) {
              B[i] = B[i - 1] + C[i] * h[i] - D[i] / 2 * h[i] * h[i];
          }
      }

      T interpolate(const T& x_0) const noexcept {
          T y_0;
          std::size_t N = x.size();

          for (std::size_t i = 0; i < N - 1; i++) {
              if (x_0 >= x[i] && x_0 <= x[i + 1]) {
                  y_0 = A[i]
                       + B[i] * (x_0 - x[i + 1])
                       + C[i] / 2 * (x_0 - x[i + 1]) * (x_0 - x[i + 1])
                       + D[i] / 6 * (x_0 - x[i + 1]) * (x_0 - x[i + 1]) * (x_0 - x[i + 1]);
                  break;
              }
          }

          return y_0;
      }
  };
}

#endif
