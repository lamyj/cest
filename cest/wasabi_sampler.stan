data {
    // Number of frequency offsets
    int<lower=0> N;
    
    // RF offset w.r.t. Larmor frequency [rad/s]
    vector[N] Delta_w;
    
    // Normalized signal magnitude
    vector[N] Z;
    
    // Main magnetic field [T]
    real<lower=0> B0;
    
    // Nominal RF amplitude [T];
    real<lower=0> B1_nominal;
    
    // RF pulse duration [s]
    real t_p;
}

transformed data {
    // Larmor constant of the proton [rad/s/T]
    // https://www.physics.nist.gov/cgi-bin/cuu/Value?gammap
    real gamma = 2.6752218708*1e8;
    
    // Main magnetic field [rad/s]
    real w0 = gamma * B0;
}

parameters {
    // Amplitude modulation independent of the frequency offset
    real<lower=0> c, d;
    // RF amplitude [T]
    real<lower=0> B1;
    // B0 inhomogeneity [rad/s]
    real delta_w;
    
    // Standard deviation of the likelihood
    real<lower=1.2e-38, upper=3.4e+38> sigma;
}

model {
    c ~ normal(0.5, 1);
    d ~ normal(1, 1);
    B1 ~ normal(B1_nominal, 1e-6);
    // TODO: How many ppms for variance?
    delta_w ~ normal(0, w0*1e-6);
    
    vector[N] w = Delta_w-delta_w;
    real w1 = gamma * B1;
    
    vector[N] mu = abs(
        c
        - d
            .* sin(atan(w1 / w)) ^ 2
            .* sin(sqrt(w1^2 + w^2) * t_p/2) ^ 2
    );
    
    sigma ~ exponential(2);
    
    Z ~ normal(mu, sigma);
}
