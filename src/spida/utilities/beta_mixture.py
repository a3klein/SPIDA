# Helper and model using component-specific mu_k and kappa priors
import jax
import jax.numpy as jnp
import numpyro.distributions as dist
import numpyro

# convert mu (in (0,1)) and kappa (>0) to alpha,beta
def alpha_beta_from_mu_kappa(mu, kappa):
    """Return (alpha, beta) for Beta distribution from mean mu and concentration kappa.
    alpha = mu * kappa, beta = (1-mu) * kappa.
    mu and kappa can be scalars or arrays (broadcastable).
    """
    mu = jnp.asarray(mu)
    kappa = jnp.asarray(kappa)
    alpha = mu * kappa
    beta = (1.0 - mu) * kappa
    return alpha, beta

# Beta mixture model variant that places priors on component-specific mu_k and kappa_k
K = 2

def beta_mixture_model_mu_kappa(scaled_data=None,
                                 mu_loc_prior=None,
                                 mu_conc_prior=(2.0, 2.0),
                                 kappa_prior_params=(2.0, 0.5)):
    """
    Marginalized Beta mixture where each component k has Beta(alpha_k, beta_k)
    with alpha_k = mu_k * kappa_k and beta_k = (1-mu_k) * kappa_k.

    Arguments:
    - scaled_data: array in (0,1)
    - mu_loc_prior: if None, we sample mu_k ~ Beta(mu_conc_prior[0], mu_conc_prior[1]) independently for each k.
                    if provided, should be an array-like of length K and will be used as fixed mu_k values.
    - mu_conc_prior: two values (a,b) used as Beta(a,b) when sampling mu_k.
    - kappa_prior_params: (shape, rate) parameters for Gamma prior on kappa_k.
    """
    weights = numpyro.sample('weights', dist.Dirichlet(0.5 * jnp.ones(K)))

    # component-specific mu_k (in (0,1))
    if mu_loc_prior is None:
        # sample each mu_k from a Beta prior with modest mass around 0.5 by default
        mu_k = numpyro.sample('mu_k', dist.Beta(jnp.ones(K) * mu_conc_prior[0], jnp.ones(K) * mu_conc_prior[1]))
    else:
        mu_k = jnp.asarray(mu_loc_prior)

    # component-specific concentration kappa_k > 0 (how tight the beta is around mu_k)
    kappa_k = numpyro.sample('kappa_k', dist.Gamma(jnp.ones(K) * kappa_prior_params[0],
                                                  jnp.ones(K) * kappa_prior_params[1]))

    # convert to alpha/beta
    alpha_k, beta_k = alpha_beta_from_mu_kappa(mu_k, kappa_k)

    # ensure data in (0,1)
    u = jnp.clip(scaled_data, 1e-6, 1.0 - 1e-6)

    # stack per-component log-probs
    comp_logp = jnp.stack([dist.Beta(alpha_k[k], beta_k[k]).log_prob(u) for k in range(K)], axis=1)

    log_w = jnp.log(weights + 1e-12)
    total_logp = jax.scipy.special.logsumexp(log_w + comp_logp, axis=1)
    numpyro.factor('obs_logp', total_logp.sum())


def to_unit(x, data_max, data_min):
    return (x - data_min) / (data_max - data_min + 1e-12)

def from_unit(u, data_max, data_min):
    return u * (data_max - data_min) + data_min

def run_beta_dist_model(data): 
    # Scale data to (0,1)
    data_min = float(np.min(np.array(data)))
    data_max = float(np.max(np.array(data)))
    scale_eps = 1e-6

    scaled = to_unit(np.array(data), data_max, data_min).astype(float)
    scaled = np.clip(scaled, scale_eps, 1.0 - scale_eps)
    scaled_jnp = jnp.array(scaled)

    # Define marginalized Beta mixture model
    K = 2

    # For my specific example These are the priors I am using! 
    mu_loc = np.percentile(scaled, [5, 50])
    kappa_params = [1.0, 5.0]  # shape, rate for Gamma prior on kappa_k

    # Run NUTS
    rng_key = random.PRNGKey(123)
    kernel = NUTS(beta_mixture_model_mu_kappa)
    mcmc_kappa = MCMC(kernel, num_warmup=800, num_samples=1200, num_chains=4)
    mcmc_kappa.run(rng_key, scaled_data=scaled_jnp, mu_loc_prior=mu_loc, kappa_prior_params=kappa_params)
    # mcmc_kappa.print_summary()

    print('alpha_beta_from_mu_kappa defined; use beta_mixture_model_mu_kappa to fit component-specific Beta priors.')

    posterior = mcmc_kappa.get_samples()
    post_mu_k = posterior.get('mu_k') # present only if mu_k was sampled
    post_kappa = posterior['kappa_k'] # shape (num_draws, K)
    post_weights = posterior['weights']

    if post_mu_k is None:
        mu_vals_for_alpha = mu_loc[np.newaxis, :] * np.ones((post_kappa.shape[0], K))
    else:
        mu_vals_for_alpha = np.array(post_mu_k)

    post_alpha = mu_vals_for_alpha * np.array(post_kappa)
    post_beta = (1.0 - mu_vals_for_alpha) * np.array(post_kappa)
    print('posterior alpha mean:', post_alpha.mean(axis=0))
    print('posterior beta mean:', post_beta.mean(axis=0))
    print("posterior post kappa: ", post_kappa.mean(axis=0))

    # Average posterior responsibilities (on original scale): use pdfs of each component
    num_resp_draws = min(1000, post_alpha.shape[0])
    sel = np.random.choice(post_alpha.shape[0], size=num_resp_draws, replace=False)
    N = len(data)
    resp_acc = np.zeros((N, K))
    for i in sel:
        w = post_weights[i]
        a = post_alpha[i]
        b = post_beta[i]
        # component pdfs on original scale
        comp0 = (1.0 / (data_max - data_min + 1e-12)) * sp_beta.pdf(to_unit(data_np, data_max, data_min), a=a[0], b=b[0])
        comp1 = (1.0 / (data_max - data_min + 1e-12)) * sp_beta.pdf(to_unit(data_np, data_max, data_min), a=a[1], b=b[1])
        weighted = np.vstack([comp0 * w[0], comp1 * w[1]]).T
        denom = weighted.sum(axis=1, keepdims=True) + 1e-12
        ratio = weighted / denom
        resp_acc += ratio

    resp_mean = resp_acc / num_resp_draws
    mod_0_support = resp_mean[:, 0]
    thr_line = data_np[(mod_0_support < 0.5).argmax()]

    return data_max, data_min, post_alpha, post_beta, post_weights


def plot_fit_beta_results(data, data_np, data_max, data_min, post_alpha, post_beta, post_weights): 



    return thr_line