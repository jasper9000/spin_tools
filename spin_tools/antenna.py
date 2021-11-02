import skrf as rf
from scipy.optimize import differential_evolution
from scipy.integrate import trapz


def optim_func_factory(line, raw_ntw, x, mask_f, limit_depth):
    def optimization_decorator(network_func):
        def wrapper(res_x):
            ntw = network_func(line, raw_ntw, res_x)
            y_f = ntw.s_db.reshape(-1)
            mask_limit_depth_f = y_f>=limit_depth

            target = trapz(y_f[mask_f&mask_limit_depth_f], x[mask_f&mask_limit_depth_f])
            target *= 1e-6
            # equal peak amplitudes
            #target += np.abs(np.min(y_f1[mask_f1]) - np.min(y_f2[mask_f2])) * 1.3
            return target
        return wrapper
    return optimization_decorator


def optimize_antenna(reduced_raw_network, component_bounds, matching_networks, f_opti_mask, limit_depth=-40, Z0=50):
    frequency_reduced = reduced_raw_network.frequency

    beta_reduced = frequency_reduced.w/rf.c
    line_reduced = rf.DefinedGammaZ0(frequency=frequency_reduced,
                                 gamma=1j*beta_reduced, z0=Z0)

    df = frequency_reduced.df[0]
    x = frequency_reduced.f

    mask_f = (x>f_opti_mask[0])&(x<f_opti_mask[1])
    results = []

    for matching_func in matching_networks:
        print(f"Optimizing {matching_func.__name__}.")
        func = optim_func_factory(line_reduced, reduced_raw_network, x, mask_f, limit_depth)(matching_func)
        res = differential_evolution(func, component_bounds,
                                    maxiter=250, popsize=30, tol=0.000001, atol=0.000001,
                                    seed=420,
                                    polish=True, disp=False)
        results.append(res)
    return results