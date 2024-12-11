import mosek
import numpy as np
import pandas as pd
import random
import time
import os
import scipy
from scipy.optimize import minimize, Bounds
import warnings
import cvxpy as cp
import math
import concurrent.futures

def approx_coeff(coeffs, state, t, x, y):
    if any(val < 0 for val in x) or y < 0:
        return -99999
    elif sum(x) == 0 or y == 0:
        return 0
    else:
        return np.dot(coeffs[t, :], state)

def getnewstate(x, y, eps, VFA, state, w, ai):
    sz = state.shape[0]
    if any(val < 0 for val in x) or y < 0:
        statenew = np.zeros(sz)
    elif np.sum(x) == 0 or y == 0:
        statenew = np.zeros(sz)
    else:
        adj = np.sum(x[0::2] + x[1::2] == 2)
        if y == 1:
            adj = 0
        if VFA == 0:
            statenew = np.array([
                1,
                np.sum(x),
                y,
                adj,
                np.sqrt(np.sum(x)),
                np.sqrt(y),
                np.sqrt(adj),
                np.sum(x)**2,
                y**2,
                adj**2])
        elif VFA == 1:
            statenew = np.array([
                1,
                np.sum(x[w]),
                np.sum(x[ai]),
                y,
                adj,
                np.sqrt(np.sum(x[w])),
                np.sqrt(np.sum(x[ai])),
                np.sqrt(y),
                np.sqrt(adj),
                np.sum(x[w]) ** 2,
                np.sum(x[ai]) ** 2,
                y ** 2,
                adj ** 2])
        elif VFA == 2:
            statenew = np.array([
                1,
                np.sum(x[:c//4]),
                np.sum(x[c//4:c//2]),
                np.sum(x[c//2:3*c//4]),
                np.sum(x[3*c//4:]),
                y,
                adj,
                np.sqrt(np.sum(x[:c//4])),
                np.sqrt(np.sum(x[c//4:c//2])),
                np.sqrt(np.sum(x[c//2:3*c//4])),
                np.sqrt(np.sum(x[3*c//4:])),
                np.sqrt(y),
                np.sqrt(adj),
                np.sum(x[:c//4]) ** 2,
                np.sum(x[c//4:c//2]) ** 2,
                np.sum(x[c//2:3*c//4]) ** 2,
                np.sum(x[3*c//4:]) ** 2,
                y ** 2,
                adj ** 2])
        elif VFA == 3:
            statenew = np.hstack([
                1,
                x,
                y,
                adj])
        elif VFA == 4:
            xstar = np.zeros(x.shape[0])
            xstar[::2] = x[1::2]
            xstar[1::2] = x[::2]
            statenew = np.hstack([
                1,
                x,
                y,
                (xstar*x)[::2]])
    return statenew

def getallstates(M,x1,Y,adj,VFA,w,ai):
    if VFA == 0:
        allstates = np.column_stack((
            np.ones(M),
            np.sum(x1, axis=0),
            Y,
            adj,
            np.sqrt(np.sum(x1, axis=0)),
            np.sqrt(Y),
            np.sqrt(adj),
            np.sum(x1, axis=0) ** 2,
            Y ** 2,
            adj ** 2))
    elif VFA == 1:
        allstates = np.column_stack([
            np.ones(M),
            np.sum(x1[w, :], axis=0),
            np.sum(x1[ai, :], axis=0),
            Y,
            adj,
            np.sqrt(np.sum(x1[w, :], axis=0)),
            np.sqrt(np.sum(x1[ai, :], axis=0)),
            np.sqrt(Y),
            np.sqrt(adj),
            np.sum(x1[w, :], axis=0) ** 2,
            np.sum(x1[ai, :], axis=0) ** 2,
            Y ** 2,
            adj ** 2])
    elif VFA == 2:
        allstates = np.column_stack((
            np.ones(M),
            np.sum(x1[:c // 4, :], axis=0),
            np.sum(x1[c // 4:c // 2, :], axis=0),
            np.sum(x1[c // 2:3 * c // 4, :], axis=0),
            np.sum(x1[3 * c // 4:, :], axis=0),
            Y,
            adj,
            np.sqrt(np.sum(x1[:c // 4, :], axis=0)),
            np.sqrt(np.sum(x1[c // 4:c // 2, :], axis=0)),
            np.sqrt(np.sum(x1[c // 2:3 * c // 4, :], axis=0)),
            np.sqrt(np.sum(x1[3 * c // 4:, :], axis=0)),
            np.sqrt(Y),
            np.sqrt(adj),
            np.sum(x1[:c // 4, :], axis=0) ** 2,
            np.sum(x1[c // 4:c // 2, :], axis=0) ** 2,
            np.sum(x1[c // 2:3 * c // 4, :], axis=0) ** 2,
            np.sum(x1[3 * c // 4:, :], axis=0) ** 2,
            Y ** 2,
            adj ** 2))
    elif VFA == 3:
        allstates = np.column_stack((
            np.ones(M),
            x1.T,
            Y,
            adj))
    elif VFA == 4:
        xstar1 = np.zeros(x1.shape)
        xstar1[::2,:] = x1[1::2,:]
        xstar1[1::2,:] = x1[::2,:]
        allstates = np.column_stack((
            np.ones(M),
            x1.T,
            Y,
            (xstar1.T*x1.T)[:,::2]))
    return allstates


def ols_algorithm_coeff(state, vmt, m):
    X = state[:m, :]

    # Exclude NaN
    valid_indices = ~np.isnan(vmt[:m])
    X = X[valid_indices, :]
    vmt2 = vmt[:m][valid_indices]
    # Compute OLS
    try:
        invers = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError as e:
        print(e)
    coeff = invers @ X.T @ vmt2

    return coeff, invers

def opt_price_BADP_coeff(c, x, y, coeffs, state, t, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3, eps, mode, VFA, w, ai):
    approxxjy1t1 = np.zeros(c)
    approxxjjy2t1 = np.zeros(c)
    x2j = x - np.eye(c)
    x3j = x - np.eye(c)

    statey1 = getnewstate(x, y - 1, eps, VFA, state, w, ai)
    approxxyt1 = approx_coeff(coeffs, state, t - 1, x, y)
    approxxy1t1 = approx_coeff(coeffs, statey1, t - 1, x, y - 1)

    for j in range(c):
        x3j[j, j - (-1) ** (j + 1)] -= 1
        approxxjy1t1[j] = approx_coeff(coeffs, getnewstate(x2j[j, :], y - 1, eps, VFA, state, w, ai), t - 1, x2j[j, :], y - 1)
        approxxjjy2t1[j] = approx_coeff(coeffs, getnewstate(x3j[j, :], y - 2, eps, VFA, state, w, ai), t - 1, x3j[j, :], y - 2)

    if mode == 'scipy':
        result = optimize_price_scipy(c, y, x, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3,
                            approxxyt1, approxxy1t1, approxxjy1t1, approxxjjy2t1)
    elif mode == 'cvx':
        result = optimize_price_cvx(c, y, x, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3,
                            approxxyt1, approxxy1t1, approxxjy1t1, approxxjjy2t1)
    ofv = result[-1] #chose last element of result vector
    return ofv

def opt_price_BADP_returnProbs_coeffs(c, x, y, coeffs, state, t, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3, eps, VFA, w, ai):
    approxxjy1t1 = np.zeros(c)
    approxxjjy2t1 = np.zeros(c)
    x2j = x - np.eye(c)
    x3j = x - np.eye(c)

    statey1 = getnewstate(x, y - 1, eps, VFA, state, w, ai)
    approxxyt1 = approx_coeff(coeffs, state, t - 1, x, y)
    approxxy1t1 = approx_coeff(coeffs, statey1, t - 1, x, y - 1)

    for j in range(c):
        x3j[j, j - (-1) ** (j + 1)] -= 1
        approxxjy1t1[j] = approx_coeff(coeffs, getnewstate(x2j[j, :], y - 1, eps, VFA, state, w, ai), t - 1, x2j[j, :], y - 1)
        approxxjjy2t1[j] = approx_coeff(coeffs, getnewstate(x3j[j, :], y - 2, eps, VFA, state, w, ai), t - 1, x3j[j, :], y - 2)

    if mode == 'scipy':
        fp1, fp2, fp3, fp2j, fp3j, r1, r2j, r3j, ofv = optimize_price_scipy(c, y, x, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3,
                            approxxyt1, approxxy1t1, approxxjy1t1, approxxjjy2t1)
    elif mode == 'cvx':
        fp1, fp2, fp3, fp2j, fp3j, r1, r2j, r3j, ofv = optimize_price_cvx(c, y, x, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3,
                            approxxyt1, approxxy1t1, approxxjy1t1, approxxjjy2t1)

    return fp1, fp2, fp3, fp2j, fp3j, r1, r2j, r3j, ofv

def optimize_price_scipy(c, y, x, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3,approxxyt1, approxxy1t1, approxxjy1t1, approxxjjy2t1):
    x = x.T
    # Bounds when all r == 0
    p0min = np.exp(a0) / (np.exp(a0) + np.exp(a1 * tau1) + np.sum(np.exp(a2) ** tau2) + np.sum(np.exp(a3) ** tau3))
    p1max = np.exp(a1 * tau1) / (np.exp(a0) + np.exp(a1 * tau1))
    p2max = np.sum(np.exp(a2) ** tau2) / (np.exp(a0) + np.sum(np.exp(a2) ** tau2))
    p3max = np.sum(np.exp(a3) ** tau3) / (np.exp(a0) + np.sum(np.exp(a3) ** tau3))

    xjstar = x.copy()
    eps = 0.00000001 * any([np.sum(x) < c - 1, np.sum(x) != y]) + 0.00000001 * all([np.sum(x) >= c - 1, np.sum(x) == y])
    xjstar[0:c:2] = x[1:c:2]
    xjstar[1:c:2] = x[0:c:2]
    adj = int(1 * (1 - (y <= 1)) * any(
        x + xjstar == 2))  # set double seat availability to 1, only if y >= 2 and double seats are available

    bounds = Bounds(np.concatenate([[eps] * 3, [eps/c] * (2 * c)]),np.concatenate(([p1max, p2max, p3max], p2max*x+eps/c, p3max*x+eps/c)))
    constraints=[]

    if adj == 1:
        def objective(vars):
            p1, p2, p3, p2j, p3j = vars[0], vars[1], vars[2], vars[3:3 + c], vars[3 + c:]

            obj = ((1.0 - p1 - p2 - p3) * approxxyt1 +
                   p1 * (a1 / b1 + 1 / (b1 * tau1) * (np.log(1 - p1 - p2 - p3) - np.log(p1)) + approxxy1t1) +
                   np.sum(p2j * (
                           a2 / b2 + 1 / b2 * (np.log(1 - p1 - p2 - p3) - np.log(p2j)) + 1 / b2 * (1 - tau2) / tau2 * (
                           np.log(1 - p1 - p2 - p3) - np.log(p2)) + approxxjy1t1)) +
                   np.sum(p3j * (
                           a3 / b3 + 1 / b3 * (np.log(1 - p1 - p2 - p3) - np.log(p3j)) + 1 / b3 * (1 - tau3) / tau3 * (
                           np.log(1 - p1 - p2 - p3) - np.log(p3)) + approxxjjy2t1)))

            return -obj  # Maximize objective by minimizing the negative

        constraints.append({'type': 'ineq', 'fun': lambda vars: x + eps / c - vars[3 + c:]})
        constraints.append({'type': 'ineq', 'fun': lambda vars: xjstar + eps / c - vars[3 + c:]})

    else:
        def objective(vars):
            p1, p2, p3, p2j, p3j = vars[0], vars[1], vars[2], vars[3:3 + c], vars[3 + c:]
            obj = ((1.0 - p1 - p2 - p3) * approxxyt1 +
                   p1 * (a1 / b1 + 1 / (b1 * tau1) * (np.log(1 - p1 - p2 - p3) - np.log(p1)) + approxxy1t1) +
                   np.sum(p2j * (
                           a2 / b2 + 1 / b2 * (np.log(1 - p1 - p2 - p3) - np.log(p2j)) + 1 / b2 * (1 - tau2) / tau2 * (
                           np.log(1 - p1 - p2 - p3) - np.log(p2)) + approxxjy1t1)) )
            return -obj  # Maximize objective by minimizing the negative

    constraints.append({'type': 'eq', 'fun': lambda vars: vars[2] - eps})
    constraints.append({'type': 'ineq', 'fun': lambda vars: 1 - p0min - vars[0] - vars[1] - vars[2]})  # nonneg
    constraints.append({'type': 'eq', 'fun': lambda vars: np.sum(vars[3:3 + c]) - vars[1]})  # sum2j
    constraints.append({'type': 'eq', 'fun': lambda vars: np.sum(vars[3 + c:]) - vars[2]})  # sum3j
    constraints.append({'type': 'ineq', 'fun': lambda vars: x + eps / c - vars[3:3 + c]})  # excludeoccupied2

    x0 = np.zeros(3 + 2 * c)
    con = (1 - p0min) / 3 - eps
    x0[0] = max(eps, con - 0.2)
    x0[1] = con + (c - np.sum(x)) * eps
    x0[3:3 + c] = (con / np.sum(x)) * x
    x0[3 + np.where(x == 0)[0]] = eps/c

    r1=47.5
    r2=47.5
    r3=47.5

    x0[0] = np.exp((a1 - b1 * r1) * tau1)/ (
                np.exp(a0) + np.exp((a1 - b1 * r1) * tau1) + sum(np.exp(a2 - b2 * r2) * tau2) + sum(
            np.exp(a3 - b3 * r3) * tau3))
    x02 = np.exp((a2 - b2* r2) * tau2) / (
                np.exp(a0) + np.exp((a1 - b1 * r1) * tau1) + sum(np.exp(a2 - b2 * r2) * tau2) + sum(
            np.exp(a3 - b3 * r3) * tau3))
    x0[3:3 + c][x == 1] = x02[x == 1]
    x0[1] = sum(x0[3:3 + c])

    if adj == 1:
        x0[2] = con + eps * (c - sum((x + xjstar == 2)))
        x0[3 + c:] = (con) / sum((x + xjstar == 2)) * (x + xjstar == 2)
        x0[3 + c:][x == 0] = eps/c
        x0[3 + c:][xjstar == 0] = eps/c
        x03 = np.exp((a3 - b3 * r3) * tau3) / (
                    np.exp(a0) + np.exp((a1 - b1 * r1) * tau1) + sum(np.exp(a2 - b2 * r2) * tau2) + sum(
                np.exp(a3 - b3 * r3) * tau3))
        x0[3 + c:][x + xjstar == 2] = x03[x + xjstar == 2]
        x0[2] = sum(x0[3 + c:])

    else:
        x0[2] = eps
        x0[3 + c:] = np.ones(c)*eps/c

    sol = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraints,options={'disp': False, 'maxiter': 10000})

    if sol.success:
        ofv = -sol.fun
        fp1, fp2, fp3 = sol.x[0], sol.x[1], sol.x[2]
        fp2j, fp3j = sol.x[3:3 + c], sol.x[3 + c:]
    else:
        ofv = float('nan')
        fp1 = fp2 = fp3 = eps
        fp2j = fp3j = np.ones(c)*(eps/c)

    r1 = a1 / b1 + 1 / (b1 * tau1) * (np.log(1.0 - fp1 - fp2 - fp3) - np.log(fp1))
    r2j = a2 / b2 + 1 / b2 * (np.log(1 - fp1 - fp2 - fp3) - np.log(fp2j)) + 1 / b2 * (1 - tau2) / tau2 * (
                np.log(1 - fp1 - fp2 - fp3) - np.log(fp2))
    r3j = a3 / b3 + 1 / b3 * (np.log(1 - fp1 - fp2 - fp3) - np.log(fp3j)) + 1 / b3 * (1 - tau3) / tau3 * (
                np.log(1 - fp1 - fp2 - fp3) - np.log(fp3))

    counter = 0
    while (r1 < -0.05 or np.any(r2j < -0.05) or np.any(r3j < -0.05)) and sol.success is False and counter < 10:
        r1 *= 0.95
        r2j *= 0.95
        r3j *= 0.95

        x0[0] = np.exp((a1 - b1 * r1) * tau1) / (
                    np.exp(a0) + np.exp((a1 - b1 * r1) * tau1) + sum(np.exp(a2 - b2 * r2) * tau2) + sum(
                np.exp(a3 - b3 * r3) * tau3))
        x02 = np.exp((a2 - b2 * r2) * tau2) / (
                    np.exp(a0) + np.exp((a1 - b1 * r1) * tau1) + sum(np.exp(a2 - b2 * r2) * tau2) + sum(
                np.exp(a3 - b3 * r3) * tau3))
        x0[3:3 + c][x == 1] = x02[x == 1]
        x0[1] = sum(x0[3:3 + c])

        if adj == 1:
            x0[2] = con + eps * (c - sum((x + xjstar == 2)))
            x0[3 + c:] = (con) / sum((x + xjstar == 2)) * (x + xjstar == 2)
            x0[3 + c:][x == 0] = eps/c
            x0[3 + c:][xjstar == 0] = eps/c
            x03 = np.exp((a3 - b3 * r3) * tau3) / (
                        np.exp(a0) + np.exp((a1 - b1 * r1) * tau1) + sum(np.exp(a2 - b2 * r2) * tau2) + sum(
                    np.exp(a3 - b3 * r3) * tau3))
            x0[3 + c:][x + xjstar == 2] = x03(x + xjstar == 2)
            x0[2] = sum(x0[3 + c:])
        else:
            x0[2] = eps
            x0[3 + c:] = np.ones(c)*eps/c

        sol = minimize(objective, x0, method='trust-exact', bounds=bounds, constraints=constraints, options={'disp': False, 'maxiter': 10000})

        if sol.success:
            ofv = -sol.fun
            fp1, fp2, fp3 = sol.x[0], sol.x[1], sol.x[2]
            fp2j, fp3j = sol.x[3:3 + c], sol.x[3 + c:]
        else:
            ofv = float('nan')
            fp1 = fp2 = fp3 = eps
            fp2j = fp3j = np.ones(c)*(eps/c)

        r1 = a1 / b1 + 1 / (b1 * tau1) * (np.log(1.0 - fp1 - fp2 - fp3) - np.log(fp1))
        r2j = a2 / b2 + 1 / b2 * (np.log(1 - fp1 - fp2 - fp3) - np.log(fp2j)) + 1 / b2 * (1 - tau2) / tau2 * (
                    np.log(1 - fp1 - fp2 - fp3) - np.log(fp2))
        r3j = a3 / b3 + 1 / b3 * (np.log(1 - fp1 - fp2 - fp3) - np.log(fp3j)) + 1 / b3 * (1 - tau3) / tau3 * (
                    np.log(1 - fp1 - fp2 - fp3) - np.log(fp3))

        counter += 1

    if (r1 < -0.05 or np.any(r2j < -0.05) or np.any(r3j < -0.05)):
        print([r1, *r2j, *r3j])

    return fp1, fp2, fp3, fp2j, fp3j, r1, r2j, r3j, ofv

def price1(p1,p2,p3,a1, a2, a3, b, tau):
    try: 
        if p1<=0:
            result=0
        else:
            result=a1/b[0]+1/b[0]*1/tau[0]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p1))
    except:
        result = 1000
    return result

def price2(j,p1,p2,p3,a1, a2, a3, b, tau):
    try: 
        if p2[j]<=0:
            result=0
        else:
            result=a2[j]/b[1]+1/b[1]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p2[j]))+1/b[1]*(1-tau[1])/tau[1]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(sum(p2)))
    except:
        result = 1000
    return result

def price3(j,p1,p2,p3,a1, a2, a3, b, tau):
    try:
        if p3[j]<=0:
            return 0
        if sum(p3)<=0:
            return 0
        result=a3[j]/b[2]+1/b[2]*(math.log(1-p1-sum(p2)
                                           -sum(p3))-math.log(p3[j]))\
               +1/b[2]*(1-tau[2])/tau[2]*(math.log(1-p1-sum(p2)
                                                   -sum(p3))-math.log(sum(p3)))
    except:
        result = 1000
    return result
    if p3[j]<=0:
        return 0
    if sum(p3)<=0:
        return 0
    result=a3[j]/b[2]+1/b[2]*(math.log(1-p1-sum(p2)
                                       -sum(p3))-math.log(p3[j]))\
           +1/b[2]*(1-tau[2])/tau[2]*(math.log(1-p1-sum(p2)
                                               -sum(p3))-math.log(sum(p3)))
    return result

def obj_func1(c,approxxyt1, approxxjjy2t1,approxxjy1t1,approxxy1t1,p0, p10, p20, p30, r1, r2j, r3j):
    obj_value = p10 * r1 + approxxy1t1 * p10 + sum(
        [p20[j] * (r2j[j] + approxxjy1t1[j]) for j in range(c)]) \
                + sum(
        [p30[j] * (r3j[j] + approxxjjy2t1[j]) for j in
         range(c)]) +p0 * approxxyt1
    return obj_value

def obj_func0(c,approxxyt1, approxxjjy2t1,approxxjy1t1,approxxy1t1,p0, p10, p20, p30, r1, r2j, r3j):
    obj_value = p10 * r1 + approxxy1t1 * p10 + sum(
        [p20[j] * (r2j[j] + approxxjy1t1[j]) for j in range(c)]) \
                + p0 * approxxyt1
    return obj_value

def optimize_price_cvx(c, y, x, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3,approxxyt1, approxxy1t1, approxxjy1t1, approxxjjy2t1):
    x = x.T
    # Bounds when all r == 0
    p0min = np.exp(a0) / (np.exp(a0) + np.exp(a1 * tau1) + np.sum(np.exp(a2) ** tau2) + np.sum(np.exp(a3) ** tau3))
    p1max = np.exp(a1 * tau1) / (np.exp(a0) + np.exp(a1 * tau1))
    p2max = np.sum(np.exp(a2) ** tau2) / (np.exp(a0) + np.sum(np.exp(a2) ** tau2))
    p3max = np.sum(np.exp(a3) ** tau3) / (np.exp(a0) + np.sum(np.exp(a3) ** tau3))

    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")

    xjstar = x.copy()
    eps = 0.00000001 * any([np.sum(x) < c - 1, np.sum(x) != y]) + 0.00000001 * all([np.sum(x) >= c - 1, np.sum(x) == y])
    xjstar[0:c:2] = x[1:c:2]
    xjstar[1:c:2] = x[0:c:2]
    adj = int(1 * (1 - (y <= 1)) * any(x + xjstar == 2))  # set double seat availability to 1, only if y >= 2 and double seats are available
    constraints=[]
    if adj == 1:
        objective_cp = 1 / (b1 * tau1) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / 1 * p1 + approxxy1t1 * p1 \
                       + cp.sum(
            [1 / b2 * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b2 * p2j[j] + approxxjy1t1[j] *
             p2j[j] for j in range(c)]) \
                       + 1 / b2 * (1 - tau2) / tau2 * (-cp.kl_div(p2, p0) - p2 + p0) \
                       + cp.sum([1 / b3 * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b3 * p3j[j] + approxxjjy2t1[j] * p3j[j] for j in range(c)]) \
                       + 1 / b3 * (1 - tau3) / tau3 * (-cp.kl_div(p3, p0) - p3 + p0) + p0 * approxxyt1
        objective = cp.Maximize(objective_cp)
        constraints=constraints+[p3j<=x+eps/c,
                                 p3j<=xjstar+eps/c]
    else:
        objective_cp = 1 / (b1 * tau1) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / 1 * p1 + approxxy1t1 * p1 \
                       + cp.sum(
            [1 / b2 * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b2 * p2j[j] + approxxjy1t1[j] *
             p2j[j] for j in range(c)]) \
                       + 1 / b2 * (1 - tau2) / tau2 * (-cp.kl_div(p2, p0) - p2 + p0) + p0 * approxxyt1
        objective = cp.Maximize(objective_cp)
        constraints=constraints+[p3==eps]

    constraints=constraints+[p1>=eps,
                             p2>=eps,
                             p3>=eps,
                             p2j>=eps/c,
                             p3j>=eps/c,
                             p1<=p1max,
                             p2<=p2max,
                             p3<=p3max,
                             p2j<=p2max*(x+eps/c),
                             p3j <= p3max * (x + eps / c),
                             p0+p1+p2+p3==1,
                             p2 == cp.sum(p2j),
                             p3 == cp.sum(p3j),
                             p0>=p0min,
                             p2j<=x+eps/c]

    subp = cp.Problem(objective, constraints)
    try:
        subp.solve(solver=cp.MOSEK)
        if math.isnan(subp.value):
            subp.solve(solver=cp.SCS)
        pass
    except cp.error.SolverError as e:
        subp.solve(solver=cp.SCS)

    fp0 = p0[0].value
    fp1 = p1[0].value
    fp2 = p2[0].value
    fp3 = p3[0].value
    fp2j = p2j.value
    fp3j = p3j.value
    ofv = subp.value
    r1 = price1(fp1, fp2j, fp3j, a1, a2, a3, [b1,b2,b3], [tau1,tau2,tau3])
    r2j = [price2(j, fp1, fp2j, fp3j, a1, a2, a3, [b1,b2,b3], [tau1, tau2,tau3]) for j in range(c)]
    r3j = [price3(j, fp1, fp2j, fp3j, a1, a2, a3, [b1,b2,b3], [tau1,tau2,tau3]) for j in range(c)]

    if ofv==float('inf') or ofv==float('-inf'):
        if adj == 1:
            ofv = obj_func1(c,approxxyt1, approxxjjy2t1, approxxjy1t1, approxxy1t1, fp0, fp1, fp2j, fp3j, r1, r2j, r3j)
        else:
            ofv = obj_func0(c,approxxyt1, approxxjjy2t1, approxxjy1t1, approxxy1t1, fp0, fp1, fp2j, fp3j, r1, r2j, r3j)

    return fp1, fp2, fp3, fp2j, fp3j, r1, r2j, r3j, ofv

def pricing_homo_seats(c):
    a0 = 0
    a1 = 0.2
    a2 = 0.4 * np.ones(c)
    a3 = 0.6 * np.ones(c)
    b1 = 0.2
    b2 = 0.4
    b3 = 0.6
    tau0 = 1
    tau1 = 1
    tau2 = 0.4
    tau3 = 0.6
    return a0, a1, a2, a3, b1, b2, b3, tau0, tau1, tau2, tau3

eps = 0
choice = 1      #0=homogeneous, 1=aisle/window, 2=heterogeneous
VFA = 0         #note that VFA 0,1,2 refer to the choice 0,1,2 respectively, other VFA might be used for all choices 
mode = 'cvx'    #alternatives: 'cvx' and 'scipy'
choice=0
c=8
a0, a1, a2, a3, b1, b2, b3, tau0, tau1, tau2, tau3 = pricing_homo_seats(c)
w = 0 
ai = 0
for a3i in range(51):
    a3=[a3i*0.1]*c
    print("a3",a3)
    approx = 0
    T = 2 * c
    nElemVFA = getallstates(1,np.ones([c,1]),c,c/2,VFA,0,0).shape[1] 
    coeffs = np.zeros((T + 1, nElemVFA))
    nsim = 100

    for M in [100]:
        tstart = time.time()
        for t in range(1, T + 1):
            # Sample pre-decision state
            print(t)
            success=False
            ran_seed=c * M * t+1
            attempt=0
            while success==False:
                np.random.seed(ran_seed)
                Y = np.random.randint(1, c+1, M)
                x1 = np.zeros((c, M))
                for m in range(M):
                    xx = np.random.permutation(c)[:np.random.randint(Y[m], c + 1)]
                    x1[xx, m] = 1

                adj = np.sum(x1[0::2, :] + x1[1::2, :] == 2, axis=0)
                adj[Y == 1] = 0
                allstates = getallstates(M,x1,Y,adj,VFA,w,ai) 

                vmt = np.zeros(M)

                for m in range(M):
                    x = x1[:, m]
                    y = Y[m]
                    ofv = opt_price_BADP_coeff(c, x, y, coeffs, allstates[m, :], t, a0, a1, a2, a3, b1, b2, b3, tau1,
                                            tau2, tau3, eps, mode, VFA, w, ai)
                    success=True
                    vmt[m] = ofv

                # Compute regression
                coeff, invers = ols_algorithm_coeff(allstates, vmt, M)
                coeffs[t, :] = coeff

        t_training = time.time() - tstart

        x = np.ones(c)
        y = np.sum(x)
        adj = np.sum(x[0::2] + x[1::2] == 2)
        adj = 0 if y == 1 else adj
        t = T
        if VFA == 0:          
            state = np.array([1, np.sum(x), y, adj, np.sqrt(np.sum(x)), np.sqrt(y), np.sqrt(adj), np.sum(x) ** 2, y ** 2, adj ** 2])
        elif VFA == 1:
            state = np.array([1,np.sum(x[w]),np.sum(x[ai]),
                y,adj,np.sqrt(np.sum(x[w])),np.sqrt(np.sum(x[ai])),
                np.sqrt(y),np.sqrt(adj),np.sum(x[w]) ** 2,
                np.sum(x[ai]) ** 2,y ** 2,adj ** 2])
        elif VFA == 2:
            state = np.array([1,np.sum(x[:c // 4]),
                np.sum(x[c // 4:c // 2]),
                np.sum(x[c // 2:3 * c // 4]),
                np.sum(x[3 * c // 4:]),y,adj,
                np.sqrt(np.sum(x[:c // 4])),
                np.sqrt(np.sum(x[c // 4:c // 2])),
                np.sqrt(np.sum(x[c // 2:3 * c // 4])),
                np.sqrt(np.sum(x[3 * c // 4:])),
                np.sqrt(y),np.sqrt(adj),
                np.sum(x[:c // 4]) ** 2,
                np.sum(x[c // 4:c // 2]) ** 2,
                np.sum(x[c // 2:3 * c // 4]) ** 2,
                np.sum(x[3 * c // 4:]) ** 2,
                y ** 2,adj ** 2])
        elif VFA == 3:
            state = np.hstack([
                1,
                x,
                y,
                adj])
        elif VFA == 4:
            xstar = np.zeros(x.shape[0])
            xstar[::2] = x[1::2]
            xstar[1::2] = x[::2]
            state = np.hstack([
                1,
                x,
                y,
                (xstar*x)[::2]])
        Vroot = approx_coeff(coeffs, state, t, x, y)

        folder_path = "EmptySeatsStudy (experiments)/sbADP"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_name = f"sbADP_learning{choice}capacity{c}a3{a3i*0.1}.mat"
        file_path = os.path.join(folder_path, file_name)
        data_to_save = {
            'coeffs': coeffs,
            't_training': t_training,
            'Vroot': Vroot,
            'a0': a0,
            'a1': a1,
            'a2': a2,
            'a3': a3,
            'b1': b1,
            'b2': b2,
            'b3': b3,
            'tau0': tau0,
            'tau1': tau1,
            'tau2': tau2,
            'tau3': tau3,
            'c': c,
            'T': T,
            'VFA': VFA,
            'w': w,
            'ai': ai}
        # Use a library like scipy.io.savemat if needed to save .mat files
        scipy.io.savemat(file_path, data_to_save)

       