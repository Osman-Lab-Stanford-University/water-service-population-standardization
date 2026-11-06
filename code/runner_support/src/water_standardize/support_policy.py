"""Predeclared diagnostic reporting screens, not identification tests."""
import math
import numpy as np

LIMITS = {
    'minimum_regime_n': 30,
    'minimum_weight_ess': 30.0,
    'minimum_ess_fraction': 0.25,
    'maximum_single_weight_fraction': 0.10,
    'minimum_target_propensity_p01': 0.001,
    'maximum_cdf_boundary_overshoot': 0.005,
    'maximum_cdf_projection_change': 0.02,
}


def _number(v, name):
    if isinstance(v, (bool, str, bytes)):
        raise ValueError(f'Invalid numeric field {name}')
    try:
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError(f'Invalid numeric field {name}') from None
    if not math.isfinite(x):
        raise ValueError(f'Nonfinite numeric field {name}')
    return x


def assess_regime(d):
    required = ('n', 'ess', 'ess_fraction', 'largest_weight_fraction', 'p01_target')
    if not isinstance(d, dict) or any(k not in d for k in required):
        raise ValueError('Incomplete overlap diagnostics')
    n = _number(d['n'], 'n')
    ess = _number(d['ess'], 'ess')
    fraction = _number(d['ess_fraction'], 'ess_fraction')
    concentration = _number(d['largest_weight_fraction'], 'largest_weight_fraction')
    p01 = _number(d['p01_target'], 'p01_target')
    if (not n.is_integer() or n < 1 or not 0 < ess <= n + 1e-6
            or not 0 < fraction <= 1 or abs(fraction - ess / n) > 1e-6
            or not 0 < concentration <= 1 or concentration < 1 / n - 1e-9
            or not 0 <= p01 <= 1):
        raise ValueError('Internally inconsistent overlap diagnostics')
    reasons = []
    if n < LIMITS['minimum_regime_n']: reasons.append('N_BELOW_30')
    if ess < LIMITS['minimum_weight_ess']: reasons.append('ESS_BELOW_30')
    if fraction < LIMITS['minimum_ess_fraction']: reasons.append('ESS_FRACTION_BELOW_0_25')
    if concentration > LIMITS['maximum_single_weight_fraction']: reasons.append('WEIGHT_CONCENTRATION_ABOVE_0_10')
    if p01 < LIMITS['minimum_target_propensity_p01']: reasons.append('TARGET_PROPENSITY_P01_BELOW_0_001')
    return {
        'screen_pass': not reasons,
        'status': 'SUPPORT_SCREEN_PASS_ONLY' if not reasons else 'SUPPORT_SCREEN_FAIL',
        'reasons': reasons,
        'identification_proven': False,
        'inference_authorized': False,
    }


def _check_screen(d):
    if (not isinstance(d, dict)
            or not isinstance(d.get('screen_pass'), bool)
            or d.get('status') != ('SUPPORT_SCREEN_PASS_ONLY' if d['screen_pass'] else 'SUPPORT_SCREEN_FAIL')
            or d.get('inference_authorized') is not False
            or d.get('identification_proven') is not False
            or not isinstance(d.get('reasons'), list)
            or bool(d['reasons']) == d['screen_pass']):
        raise ValueError('Malformed or bypassed group screen')


def assess_pair(a, b):
    _check_screen(a)
    _check_screen(b)
    permitted = a['screen_pass'] and b['screen_pass']
    return {
        'screen_pass': permitted,
        'status': 'PAIR_SUPPORT_SCREEN_PASS_ONLY' if permitted else 'PAIR_SUPPORT_SCREEN_FAIL',
        'inference_authorized': False,
    }


def assess_curve(raw, projected, grid, support):
    _check_screen(support)
    raw = np.asarray(raw, dtype=float)
    projected = np.asarray(projected, dtype=float)
    grid = np.asarray(grid, dtype=float)
    if (raw.ndim != 1 or projected.shape != raw.shape or grid.shape != raw.shape
            or len(raw) < 2 or not np.isfinite(raw).all()
            or not np.isfinite(projected).all() or not np.isfinite(grid).all()
            or grid[0] != 0 or grid[-1] != 100 or np.any(np.diff(grid) <= 0)
            or abs(raw[-1] - 1.0) > 1e-10
            or abs(projected[-1] - 1.0) > 1e-10
            or np.any(projected < 0) or np.any(projected > 1)
            or np.any(np.diff(projected) < -1e-10)):
        raise ValueError('Malformed CDF or endpoint violation')
    overshoot = float(max(np.maximum(-raw, 0).max(), np.maximum(raw - 1, 0).max()))
    change = float(np.max(np.abs(raw - projected)))
    steps = int(np.sum(np.diff(raw) < -1e-10))
    reasons = []
    if not support['screen_pass']: reasons.append('REGIME_SUPPORT_FAILED')
    if overshoot > LIMITS['maximum_cdf_boundary_overshoot']: reasons.append('CDF_BOUNDARY_OVERSHOOT')
    if change > LIMITS['maximum_cdf_projection_change']: reasons.append('CDF_PROJECTION_EXCESSIVE')
    return {
        'screen_pass': not reasons,
        'status': 'CDF_SCREEN_PASS_SHAPE_SENSITIVITY_ONLY' if not reasons else 'CDF_SCREEN_FAIL',
        'reasons': reasons,
        'raw_boundary_overshoot': overshoot,
        'raw_monotonicity_steps': steps,
        'max_projection_change': change,
        'raw_preserved': True,
        'simultaneous_inference_authorized': False,
    }


def assess_risk_at_threshold(t, value, se, support):
    if t not in (50, 60, 70):
        raise ValueError('Threshold not prespecified')
    _check_screen(support)
    estimate = _number(value, 'risk_estimate')
    error = _number(se, 'risk_standard_error')
    if error < 0:
        raise ValueError('Negative standard error')
    reasons = []
    if not support['screen_pass']: reasons.append('REGIME_SUPPORT_FAILED')
    if not 0 <= estimate <= 1: reasons.append('RISK_ESTIMATE_OUTSIDE_0_1')
    return {
        'screen_pass': not reasons,
        'status': 'RISK_POINT_SCREEN_PASS_ONLY' if not reasons else 'RISK_POINT_SCREEN_FAIL',
        'reasons': reasons,
        'inference_authorized': False,
    }
