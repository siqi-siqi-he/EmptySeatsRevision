import scipy.io
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
import time

def load_mat_file(file_name):
    return scipy.io.loadmat(file_name)

def save_mat_file(file_name, data):
    scipy.io.savemat(file_name, data)

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

def simulation_coeffs_V5(nsim, T, c, coeffs, state, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3, eps, VFA, w, ai):
    rev_nt = np.zeros((nsim, T))
    policy = np.zeros((nsim, T, 2 * c + 1))
    tickets_n = np.zeros((nsim, 4))
    seats3_n = np.zeros((nsim, c))
    seats2_n = np.zeros((nsim, c))

    for n in range(nsim):
        x = np.ones(c)
        y = np.sum(x)
        cTickets = np.zeros(4)
        cSeats2 = np.zeros(c)
        cSeats3 = np.zeros(c)
        adj = np.sum(x[0::2] + x[1::2] == 2)
        if y == 1:
            adj = 0
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
                
        rev_nt_local = np.zeros(T)
        policy_local = np.zeros((T, 2 * c + 1))

        for t in range(T):
            t2 = T + 1 - t
            if np.sum(x) == 0 or y == 0:
                break
            else:
                fp1, fp2, fp3, fp2j, fp3j, r1, r2j, r3j, ofv = opt_price_BADP_returnProbs_coeffs(
                    c, x, y, coeffs, state, t2, a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3, eps, VFA, w, ai)
                policy_local[t, :] = np.concatenate(([r1], r2j, r3j))
                x, y, rev, cTickets, cSeats2, cSeats3 = calc_decision_V5(fp1, fp2j, fp3j, r1, r2j, r3j, x, y, t, n, cTickets, cSeats2, cSeats3)
                state = getnewstate(x, y, eps, VFA, state, w, ai)
                rev_nt_local[t] = max(0, rev)

        rev_nt[n] = rev_nt_local
        policy[n] = policy_local
        tickets_n[n] = cTickets
        seats3_n[n] = cSeats3
        seats2_n[n] = cSeats2
        #print("revenue at simulation: ", n, ": ", rev_nt[n])

    return rev_nt, policy, tickets_n, seats2_n, seats3_n

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


def calc_decision_V5(fp1, fp2j, fp3j, r1, r2j, r3j, x, y, t, n, cTickets, cSeats2, cSeats3):
    folder_path = "EmptySeatsStudy (experiments)/RandomNumbers"
    file_name = "randomdecisions_capacity8_choice0_sim100_a3_4_t" + str(t+1) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    dec = np.loadtxt(file_path)
    # Let the customer decide
    if dec[n] > (fp1 + np.sum(fp2j) + np.sum(fp3j)):
        # Leave without a ticket
        cTickets[0] += 1
        rev = 0
        return x, y, rev, cTickets, cSeats2, cSeats3
    else:
        concat = np.concatenate(([fp1], fp2j.flatten(), fp3j.flatten()))
        cumsumConcat = np.cumsum(concat)
        lower = dec[n] < cumsumConcat
        ind = np.where(lower)[0][0]

        if ind == 0:
            # Buy product 1 (without reservation)
            y -= 1
            rev = r1
            cTickets[1] += 1
        elif ind <= len(fp2j) :
            # Buy product 2 at seat ind
            x[ind-1 ] = 0
            rev = r2j[ind-1 ]
            y -= 1
            cTickets[2] += 1
            cSeats2[ind -1] += 1
        else:
            # Buy product 3 at seat ind and its adjacent one
            ind -= len(fp2j) + 1
            x[ind] = 0
            adj_ind = ind - (-1) ** (ind+1)
            x[adj_ind] = 0
            rev = r3j[ind]
            y -= 2
            cSeats3[ind] += 1
            cSeats3[adj_ind] += 1
            cTickets[3] += 1

    return x, y, rev, cTickets, cSeats2, cSeats3

def sample_dec(T, *args):
    sz = [T]
    for arg in args:
        sz.append(arg)
    dec = np.random.rand(*sz)
    return dec

def find_empty_row(file_path, sheet):
    df = pd.read_excel(file_path, sheet)
    index = df.index[df.isnull().all(axis=1)].tolist()
    if index:
        return index[0] + 2
    else:
        return len(df) + 1

# Load data
M = 100
nsim = 100
choice = 2
VFA = 0
c=8
mode = 'cvx'
choice=0
M=100
rev_a3j=[0]*51
stderr_a3j=[0]*51
seats3_a3j=[0]*51
for a3j in range(51):
    a3=[a3j*0.1]*c
    folder_path = "EmptySeatsStudy (experiments)/sbADP"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = f"sbADP_learning{choice}capacity{c}a3{a3j*0.1}.mat"
    file_path = os.path.join(folder_path, file_name)
    mat_data = load_mat_file(file_path)

    coeffs = mat_data['coeffs']
    t_training = mat_data['t_training']
    Vroot = mat_data['Vroot']
    a0 = mat_data['a0'][0][0]
    a1 = mat_data['a1'][0][0]
    a2 = mat_data['a2'][0]
    # a3 = mat_data['a3'][0]
    b1 = mat_data['b1'][0][0]
    b2 = mat_data['b2'][0][0]
    b3 = mat_data['b3'][0][0]
    tau0 = mat_data['tau0'][0][0]
    tau1 = mat_data['tau1'][0][0]
    tau2 = mat_data['tau2'][0][0]
    tau3 = mat_data['tau3'][0][0]
    T = mat_data['T'][0][0]  # Extract scalar from nested array
    eps = 0  # Extract scalar from nested array
    w = mat_data['w'][0][0]
    ai = mat_data['ai'][0][0]

    rev_nt, policy, tickets_n, seats2_n, seats3_n = simulation_coeffs_V5(nsim, T, c, coeffs, [], a0, a1, a2, a3, b1, b2, b3, tau1, tau2, tau3, eps, VFA, w, ai)
    # tend = time.time() - tstart
    print(rev_nt)
    rev_a3j[a3j] = np.mean(np.sum(rev_nt, axis=1))
    stderr_a3j[a3j] = np.std(np.sum(rev_nt, axis=1))
    seats3_a3j[a3j] = np.mean(np.sum(seats3_n, axis=1))

folder_path = "EmptySeatsStudy (experiments)/sbADP"
file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, rev_a3j)
file_name = "simu_vars_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, stderr_a3j)
file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, seats3_a3j)
