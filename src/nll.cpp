#include <TMB.hpp>
#include <cmath>

template<class Type>
Type two_pi() { return Type(2.0) * M_PI; }

// von Mises density on [0, 2*pi). Uses scaled Bessel for stability
template<class Type>
Type von_mises_pdf(Type x, Type mu, Type kappa) {
  // density = exp(kappa * cos(x-mu)) / (2*pi*I0(kappa))
  Type I0 = besselI(kappa, Type(0));
  return exp(kappa * CppAD::cos(x - mu)) / (two_pi<Type>() * I0);
}

// Deterministic mapping from true age (years) to expected annuli count
template<class Type>
int calculate_K_for_age(Type age, Type survey_date, Type annuli_date, Type annuli_min_age) {
  if (age < Type(0)) return 0;
  Type birth_date = survey_date - age;
  Type birth_year = Type(std::floor(asDouble(birth_date)));
  Type next_ring_date = birth_year + annuli_date;
  if (next_ring_date <= birth_date) {
    next_ring_date = (birth_year + Type(1)) + annuli_date;
  }
  int k_count = 0;
  while (next_ring_date < survey_date) {
    Type age_at_ring_date = next_ring_date - birth_date;
    if (age_at_ring_date >= annuli_min_age) {
      k_count += 1;
    }
    next_ring_date = next_ring_date + Type(1);
  }
  return k_count;
}

template<class Type>
Type objective_function<Type>::operator() () {
  // Observation data
  DATA_IVECTOR(obs_survey_index);
  DATA_IVECTOR(obs_length_index);
  DATA_IVECTOR(obs_K);
  DATA_VECTOR(obs_count);

  // Survey meta
  DATA_VECTOR(survey_dates);

  // Grids
  DATA_VECTOR(l_grid);
  DATA_VECTOR(a_grid);

  // Spawning and observation constants
  DATA_SCALAR(spawning_mu);
  DATA_SCALAR(spawning_kappa);
  DATA_SCALAR(annuli_date);
  DATA_SCALAR(annuli_min_age);

  // Numerical constants
  DATA_SCALAR(Delta_l);
  DATA_SCALAR(Delta_t);
  DATA_SCALAR(log_eps);

  // Parameters
  PARAMETER(k);
  PARAMETER(L_inf);
  PARAMETER(d);
  PARAMETER(m);

  int N_l = l_grid.size();
  int N_t = a_grid.size() - 1;

  // Interfaces
  vector<Type> l_interfaces(N_l + 1);
  for (int i = 0; i <= N_l; ++i) l_interfaces(i) = Type(i) * Delta_l;

  // Coefficients
  vector<Type> v(N_l + 1), D(N_l + 1);
  for (int i = 0; i <= N_l; ++i) {
    v(i) = k * (L_inf - l_interfaces(i)) - d / Type(2.0);
    D(i) = d * l_interfaces(i) / Type(2.0);
  }
  vector<Type> v_plus(N_l + 1), v_minus(N_l + 1);
  for (int i = 0; i <= N_l; ++i) {
    v_plus(i)  = CppAD::CondExpGt(v(i), Type(0), v(i), Type(0));
    v_minus(i) = CppAD::CondExpLt(v(i), Type(0), v(i), Type(0));
  }
  vector<Type> mu_vec(N_l);
  for (int i = 0; i < N_l; ++i) mu_vec(i) = m / l_grid(i);

  // Tridiagonal system
  vector<Type> a_(N_l - 1), b_(N_l), c_(N_l - 1);
  Type c1 = Delta_t / Delta_l;
  Type c2 = Delta_t / (Delta_l * Delta_l);
  for (int i = 1; i <= N_l - 2; ++i) {
    a_(i - 1) = -c1 * v_plus(i) - c2 * D(i);
    c_(i)     =  c1 * v_minus(i + 1) - c2 * D(i + 1);
    b_(i)     =  Type(1) + Delta_t * mu_vec(i)
               + c1 * (v_plus(i + 1) - v_minus(i))
               + c2 * (D(i + 1) + D(i));
  }
  b_(0) = Type(1) + Delta_t * mu_vec(0) + c1 * (v_plus(1) + D(1) / Delta_l);
  c_(0) = c1 * (v_minus(1) - D(1) / Delta_l);
  a_(N_l - 2) = -c1 * v_plus(N_l - 1) - c2 * D(N_l - 1);
  b_(N_l - 1) = Type(1) + Delta_t * mu_vec(N_l - 1)
              + c1 * (v_plus(N_l) - v_minus(N_l - 1))
              + c2 * (D(N_l) + D(N_l - 1));

  // Thomas solver (AD friendly)
  auto solve_tridiag = [&](const vector<Type>& a, const vector<Type>& b,
                           const vector<Type>& c, const vector<Type>& dvec) {
    int n = b.size();
    vector<Type> c_prime(n);
    vector<Type> d_prime(n);
    c_prime(0) = c(0) / b(0);
    d_prime(0) = dvec(0) / b(0);
    for (int i = 1; i <= n - 2; ++i) {
      Type mloc = b(i) - a(i - 1) * c_prime(i - 1);
      c_prime(i) = c(i) / mloc;
      d_prime(i) = (dvec(i) - a(i - 1) * d_prime(i - 1)) / mloc;
    }
    d_prime(n - 1) = (dvec(n - 1) - a(n - 2) * d_prime(n - 2)) /
                     (b(n - 1) - a(n - 2) * c_prime(n - 2));
    vector<Type> x(n);
    x(n - 1) = d_prime(n - 1);
    for (int i = n - 2; i >= 0; --i) x(i) = d_prime(i) - c_prime(i) * x(i + 1);
    return x;
  };

  // Green's function matrix G: rows time (0..N_t), cols size (N_l)
  matrix<Type> G(N_t + 1, N_l);
  for (int j = 0; j < N_l; ++j) G(0, j) = Type(0);
  G(0, 0) = Type(1);
  for (int n = 0; n < N_t; ++n) {
    vector<Type> rhs(N_l);
    for (int j = 0; j < N_l; ++j) rhs(j) = G(n, j);
    vector<Type> next = solve_tridiag(a_, b_, c_, rhs);
    for (int j = 0; j < N_l; ++j) G(n + 1, j) = next(j);
  }

  int nSurvey = survey_dates.size();
  int nObs = obs_count.size();

  // Spawning weights per survey and age
  matrix<Type> spawn_w(N_t + 1, nSurvey);
  Type mu_rad = spawning_mu * two_pi<Type>();
  for (int s = 0; s < nSurvey; ++s) {
    for (int n = 0; n <= N_t; ++n) {
      Type birth_date = survey_dates(s) - a_grid(n);
      Type day_fraction = birth_date - Type(std::floor(asDouble(birth_date)));
      Type day_rad = day_fraction * two_pi<Type>();
      spawn_w(n, s) = von_mises_pdf(day_rad, mu_rad, spawning_kappa);
    }
  }

  // K(a) per survey (depends on annuli_min_age)
  matrix<int> K_of_age(N_t + 1, nSurvey);
  for (int s = 0; s < nSurvey; ++s) {
    for (int n = 0; n <= N_t; ++n) {
      K_of_age(n, s) = calculate_K_for_age(a_grid(n), survey_dates(s), annuli_date, annuli_min_age);
    }
  }

  // Denominator per (survey, length): sum over ages spawn_w * G
  matrix<Type> denom_sum(nSurvey, N_l);
  for (int s = 0; s < nSurvey; ++s) {
    for (int j = 0; j < N_l; ++j) denom_sum(s, j) = Type(0);
    for (int n = 0; n <= N_t; ++n) {
      Type w = spawn_w(n, s);
      for (int j = 0; j < N_l; ++j) denom_sum(s, j) += w * G(n, j);
    }
  }

  // NLL accumulation over observed cells only
  Type nll = Type(0);
  for (int idx = 0; idx < nObs; ++idx) {
    int s = obs_survey_index(idx) - 1;
    int j = obs_length_index(idx) - 1;
    int Kobs = obs_K(idx);
    Type num = Type(0);
    for (int n = 0; n <= N_t; ++n) if (K_of_age(n, s) == Kobs) num += spawn_w(n, s) * G(n, j);
    Type den = denom_sum(s, j);
    Type p = num / den;
    Type logp = CppAD::log(p + CppAD::exp(log_eps));
    nll -= obs_count(idx) * logp;
  }

  ADREPORT(k);
  ADREPORT(L_inf);
  ADREPORT(d);
  ADREPORT(m);

  return nll;
}


