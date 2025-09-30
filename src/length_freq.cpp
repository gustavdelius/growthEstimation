#include <TMB.hpp>
#include <cmath>

// Helper function: Sigmoidal selectivity
template<class Type>
vector<Type> calculate_sigmoid_selectivity(Type l50, Type ratio,
                              vector<Type> l)
{
    Type c1 = Type(1.0);
    Type sr = l50 * (c1 - ratio);
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;

    vector<Type> Selectivity = 1/(c1 + exp(s1 - s2 * l));

    // Ensure all elements are finite and >= 0
    TMBAD_ASSERT((Selectivity.array().isFinite() && (Selectivity.array() >= 0)).all());

    return Selectivity;
}

// Helper function: Knife edge selectivity
template<class Type>
vector<Type> calculate_knife_edge_selectivity(Type l50, Type ratio,
                                         vector<Type> l)
{
    Type c1 = Type(1.0);
    Type sr = l50 * (1.0 - 0.95); // Compute a very steep slope using ~95% of l50 as reference
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;

    vector<Type> Selectivity = 1 / (c1 + exp(s1 - s2 * l));

    // Ensure all elements are finite and >= 0
    TMBAD_ASSERT((Selectivity.array().isFinite() && (Selectivity.array() >= 0)).all());

    return Selectivity;
}

template<class Type>
Type two_pi() { return Type(2.0) * M_PI; }

// von Mises density on [0, 2*pi). Uses scaled Bessel for stability
template<class Type>
Type von_mises_pdf(Type x, Type mu, Type kappa) {
    // density = exp(kappa * cos(x-mu)) / (2*pi*I0(kappa))
    Type I0 = besselI(kappa, Type(0));
    return exp(kappa * CppAD::cos(x - mu)) / (two_pi<Type>() * I0);
}

template<class Type>
Type objective_function<Type>::operator() () {
    // Observation data
    DATA_MATRIX(count_matrix);

    // Survey meta
    DATA_VECTOR(survey_dates);

    // Grids
    DATA_VECTOR(l_grid);
    DATA_VECTOR(a_grid);

    // Spawning and observation constants
    DATA_SCALAR(spawning_mu);
    DATA_SCALAR(spawning_kappa);

    // Numerical constants
    DATA_SCALAR(Delta_l);
    DATA_SCALAR(Delta_t);
    DATA_SCALAR(log_eps);

    //Gear selectivity
    DATA_INTEGER(gear_type); // 0 = sigmoidal, 1 = knife-edge

    //selectivty params
    PARAMETER(l50);
    PARAMETER(ratio); // only used for sigmoidal gear

    //life history parameters
    PARAMETER(k);
    PARAMETER(L_inf);
    PARAMETER(d);
    PARAMETER(m);

    int N_l = l_grid.size();
    int N_t = a_grid.size() - 1;

    //Calculate number density

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

    //Calculate selectivity vector
    vector<Type> Selectivity;
    if(gear_type == 0){
        Selectivity = calculate_sigmoid_selectivity(l50, ratio, l_grid);
    } else {
        Selectivity = calculate_knife_edge_selectivity(l50, ratio, l_grid);
    }

    // Initialize negative log-likelihood accumulator
    Type nll = 0.0;

    int nSurvey = survey_dates.size();

    // Build age grid (time since birth)
    vector<Type> a_grid_shifted(N_t + 1);
    for (int i = 0; i <= N_t; ++i) {
        a_grid_shifted(i) = a_grid(i);
    }

    // Loop over each survey date
    for (int s = 0; s < nSurvey; ++s) {
        Type survey_date = survey_dates(s);

        //Calculate spawning weights for ages at survey s
        vector<Type> birth_dates(N_t + 1);
        for (int n = 0; n <= N_t; ++n) {
            birth_dates(n) = survey_date - a_grid_shifted(n);
        }

        vector<Type> spawning_weights(N_t + 1);
        Type mu_rad = spawning_mu * two_pi<Type>();
        for (int n = 0; n <= N_t; ++n) {
            Type day_fraction = birth_dates(n) - floor(asDouble(birth_dates(n)));
            Type day_rad = day_fraction * two_pi<Type>();
            spawning_weights(n) = von_mises_pdf(day_rad, mu_rad, spawning_kappa);
        }

        // Calculate population at survey s: convolve spawning weights and G
        // P(l, t) = sum all as together->N_t {G(a, l)* S(s - a)}
        // Calculate population at survey s: convolve spawning weights and G
        // G: (N_t+1) x N_l  (rows = ages/time, cols = length bins)
        // N_pop: column vector length N_l
        vector<Type> N_pop(N_l);
        N_pop.setZero();

        // Use .row(n).transpose() to produce a column vector, and
        // use .array() to ensure element-wise addition (no array/matrix mix).
        for (int n = 0; n <= N_t; ++n) {
            N_pop.array() += (spawning_weights(n) * G.row(n).transpose().array());
        }
        N_pop *= Delta_t;

        // Apply selectivity elementwise (use .array() or cwiseProduct)
        vector<Type> V_pop = N_pop.array() * Selectivity.array();

        // Normalize to get probabilities
        Type sum_V = V_pop.sum() + Type(1e-12);
        vector<Type> prob = V_pop / sum_V;

        // Calculate multinomial negative log-likelihood for survey s
        vector<Type> counts(N_l);
        for (int l = 0; l < N_l; ++l) {
            counts(l) = count_matrix(s, l);
        }
        nll -= dmultinom(counts, prob, true);
    }

    ADREPORT(k);
    ADREPORT(L_inf);
    ADREPORT(d);
    ADREPORT(m);

    // Return total negative log-likelihood summed over all surveys
    return nll;
}
