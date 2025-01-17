import matplotlib.pyplot as plt
import docplex
from docplex.mp.model import Model
import cases as cases
import math
import time
import numpy as np
import cvxpy as cp
import multiprocessing as multip
import warnings
# import mosek
import copy
from heapq import *
import itertools
import sys
from multiprocessing import Pool, Manager
import os
import pandas as pd

counter = itertools.count()
warnings.simplefilter("error", RuntimeWarning)
eps = 1e-5

class Branch_Bound():
    #for part 3
    def __init__(self, constraints=[], objective=0, x_vars=[], p1_vars=[], p2j_vars=[], p3j_vars=[],
                 y_vars=[], bar=-100):
        self.constraints = constraints
        self.objective = objective
        self.x_vars = x_vars
        self.y_vars = y_vars
        self.p1_vars = p1_vars
        self.p2j_vars = p2j_vars
        self.p3j_vars = p3j_vars
        self.children = []
        self.res_r = 1e20
        self.bar=bar

    def buildProblem(self):
        prob = cp.Problem(cp.Maximize(self.objective),
                          self.constraints)
        return prob

    def is_integral(self):
        xall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.x_vars])
        yall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.y_vars])
        return xall * yall

    def compute_obj(self, x_vars_r, y_vars_r):
        temp=copy.deepcopy(self)
        temp_x=self.x_vars
        temp_y=self.y_vars
        for j in range(len(x_vars_r)):
            temp.constraints.append(temp_x[j]==x_vars_r[j])
        temp.constraints.append(temp_y == y_vars_r[0])
        prob_temp = cp.Problem(cp.Maximize(temp.objective),
                          temp.constraints)
        try:
            res = prob_temp.solve(solver=cp.MOSEK,verbose=False)
        except cp.error.SolverError as e:
            print(e)
            try:
                res = prob_temp.solve(solver=cp.SCS, verbose=False)
            except cp.error.SolverError as e:
                print(e)
                res = prob_temp.solve(solver=cp.ECOS, verbose=False)

        return res, prob_temp

    def branch(self):
        children = []
        for branch in [0, 1]:
            child = copy.deepcopy(self)
            xmin = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.x_vars)])[0]
            ymin = [abs(v.value - 0.5) for v in self.y_vars][0]
            if abs(ymin-0.5)<1e-4 and abs(xmin-0.5)<1e-4:
                return children
            if xmin < ymin:
                r = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.x_vars)])[2]
                child.constraints.append(r == branch)
            else:
                r = self.y_vars
                child.constraints.append(r == branch)

            children.append(child)
        return children

    def bbsolve(self,time_limit=False,limit=60):
        root = self
        try:
            res = root.buildProblem().solve(solver=cp.MOSEK,verbose=False)
        except cp.error.SolverError as e:
            print(e)
            try:
                res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
            except cp.error.SolverError as e:
                print(e)
                res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)

        heap = [(-res, next(counter), root)]
        bestres_r = 1e20
        nodecount = 0
        start=time.time()
        while len(heap) > 0:
            nodecount += 1
            _, _, node = heappop(heap)
            prob = node.buildProblem()
            try:
                res = prob.solve(solver=cp.MOSEK,verbose=False)
            except cp.error.SolverError as e:
                print(e)
                try:
                    res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
                except cp.error.SolverError as e:
                    print(e)
                    res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)
            if prob.status in ["infeasible","unbounded"]:
                print("Nodes searched: ", nodecount, "infeasible")
                return bestres_r, p10, p20, p30, x0, y0

            if node.is_integral():
                if res>=bestres_r:
                    bestres_r = node.res_r
                    p10 = [v.value for v in prob.var_dict['p1']]
                    p20 = [v.value for v in prob.var_dict['p2j']]
                    p30 = [v.value for v in prob.var_dict['p3j']]
                    x0 = [v.value for v in prob.var_dict['x']]
                    y0 = [v.value for v in prob.var_dict['y']]
                continue

            if bestres_r < 1e15 and res<bestres_r:
                continue
            y_vars_r = [np.round(v.value) for v in node.y_vars]
            x_vars_r=np.array([np.round(v.value) for v in node.x_vars])
            x_vars_r[np.argpartition(np.abs(x_vars_r - 1), int(y_vars_r[0]) + 1)[:int(y_vars_r[0]) + 1]] = 1
            temp = copy.deepcopy(node)
            node.res_r, node_temp = temp.compute_obj(x_vars_r, y_vars_r)

            if bestres_r > 1e15 or bestres_r<node.res_r:
                bestres_r = node.res_r
                p10 = [v.value for v in prob.var_dict['p1']]
                p20 = [v.value for v in prob.var_dict['p2j']]
                p30 = [v.value for v in prob.var_dict['p3j']]
                x0 = x_vars_r
                y0 = y_vars_r
            if res<=self.bar:
                continue
            new_nodes = node.branch()
            for new_node in new_nodes:
                heappush(heap, (-res, next(counter),
                                new_node))
            if time_limit==True and time.time()-start>=limit:
                print("Nodes searched: ", nodecount)
                print("time:", time.time() - start)
                print("time exceeded")
                return bestres_r, p10, p20, p30, x0, y0
        end=time.time()
        print("Nodes searched: ", nodecount)
        print("time:", end-start)
        print("normal")
        return bestres_r, p10, p20, p30, x0, y0


    def solve_exist(self):
        prob = self
        try:
            res = prob.buildProblem().solve(solver=cp.MOSEK,verbose=False)
        except cp.error.SolverError as e:
            print(e)
            try:
                res = prob.buildProblem().solve(solver=cp.SCS, verbose=False)
            except cp.error.SolverError as e:
                print(e)
                res = prob.buildProblem().solve(solver=cp.ECOS, verbose=False)
        y_vars_r = [np.round(v.value) for v in prob.y_vars]
        x_vars_r = np.array([np.round(v.value) for v in prob.x_vars])
        x_vars_r[np.argpartition(np.abs(x_vars_r - 1), int(y_vars_r[0])+1)[:int(y_vars_r[0])+1]]=1
        res_r, node_temp = prob.compute_obj(x_vars_r, y_vars_r)
        p10 = [v.value for v in node_temp.var_dict['p1']]
        p20 = [v.value for v in node_temp.var_dict['p2j']]
        p30 = [v.value for v in node_temp.var_dict['p3j']]
        x0 = x_vars_r
        y0 = y_vars_r
        return res_r, p10, p20, p30, x0, y0


def UP_t(c, a1, a2, a3, b, tau):
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def UP_1(c, a1, a2, a3, b, tau):
    result = math.exp(tau[0] * (a1)) / (math.exp(tau[0] * (a1)) + 1)
    return result

def kl_div(x, y):
    return x * (math.log(x) - math.log(y)) - x + y

def UP_2(j, c, a1, a2, a3, b, tau):
    result = math.exp(a2[j]) * sum(math.exp(a2[j]) for j in range(c)) ** (tau[1] - 1) / (
                1 + sum(math.exp(a2[j]) for j in range(c)) ** (tau[1]))
    return result

def UP_3(j, c, a1, a2, a3, b, tau):
    result = math.exp(a3[j]) * sum(math.exp(a3[j]) for j in range(c)) ** (tau[2] - 1) / (
            1 + sum(math.exp(a3[j]) for j in range(c)) ** (tau[2]))
    return result

def r1(p1, p2, p3):
    if p1 <= 0:
        result = 0
    else:

        result = a1 / b[0] + 1 / b[0] * 1 / tau[0] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p1))
    return result

def r2(j, p1, p2, p3):
    if p2[j] <= 0:
        result = 0
    elif sum(p2) <= 0:
        for i in range(len(p2)):
            if p2[i] <= 0:
                p2[i] = 0
        result = a2[j] / b[1] + 1 / b[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p2[j])) + 1 / b[1] * (
                1 - tau[1]) / tau[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p2)))
    else:
        result = a2[j] / b[1] + 1 / b[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p2[j])) + 1 / b[1] * (
                    1 - tau[1]) / tau[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p2)))

    return result

def r3(j, p1, p2, p3):
    if p3[j] <= 0:
        result = 0
    elif sum(p3) <= 0:
        for i in range(len(p3)):
            if p3[i] <= 0:
                p3[i] = 0
        result = a3[j] / b[2] + 1 / b[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p3[j])) + 1 / b[2] * (
                    1 - tau[2]) / tau[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p3)))
    else:

        result = a3[j] / b[2] + 1 / b[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p3[j])) + 1 / b[2] * (
                    1 - tau[2]) / tau[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p3)))

    return result

def obj_func(v, w, theta, t, p10, p20, p30, x, y):
    obj_value = p10 * r1(p10, p20, p30) - w[t - 1] * p10 + sum(
        [p20[j] * (r2(j, p10, p20, p30) - v[t - 1, j] - w[t - 1]) for j in range(c)]) \
                + sum(
        [p30[j] * (r3(j, p10, p20, p30) - v[t - 1, j] - v[t - 1, j - (-1) ** (j + 1)] - 2 * w[t - 1]) for j in
         range(c)]) \
                - theta[t] + theta[t - 1] - (w[t] - w[t - 1]) * (y + 1) - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
    return obj_value

def subproblem(v, w, theta, t,c, b, tau, a1,a2, a3,till_end=False,time_limit=True, limit=5):
    # declare/define the decision variables of the subproblem
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    x = cp.Variable(c, name="x")
    y = cp.Variable(1, name="y")

    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 - w[t - 1] * p1 \
                   + cp.sum(
        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (v[t - 1, j] + w[t - 1]) * p2j[j]
         for j in range(c)]) \
                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                   - theta[t] + theta[t - 1] - (w[t] - w[t - 1]) * (y + 1) - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))

    constraints = [p1 <= 1,
                   p1 >= eps,
                   p0 >= eps,
                   p2 >= eps,
                   p2 <= 1,
                   p1 + p0 + p2 + p3 == 1,
                   p2 == cp.sum(p2j),
                   p3 == cp.sum(p3j),
                   p0 >= 1 - UP_t(c, a1, a2, a3, b, tau),
                   y <= 1,
                   y >= 0,
                   x >= 0,
                   x <= 1,
                   p1 <= y + 1,
                   p2j <= y + 1,
                   p3j == eps / c,  # p3j <= math.floor(y / 2) + eps / c,
                   p2j <= 1,
                   p2j >= eps / c,
                   p2j <= x + eps / c,
                   cp.sum(x)>=y+1]
    for j in range(c):
        constraints = constraints + [p3j[j] <= x[j - (-1) ** (j + 1)] + eps / c]
    bar = -100
    root = Branch_Bound(constraints=constraints, objective=objective_cp, x_vars=[x[j] for j in range(c)], p1_vars=p1,
                        p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)},
                        y_vars=y, bar=bar)
    if till_end == True:
        res, p10, p20, p30, x0, y0 = root.bbsolve(time_limit, limit)
    elif till_end == False:
        res, p10, p20, p30, x0, y0= root.solve_exist()

    return p10[0], p20, p30, x0, y0[0]+1, res

def save(v, w, theta, c, choice):
    folder_path = "wo3_BB"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, v)
    file_name = "w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, w)
    file_name = "theta_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, theta)


def execute_subproblem(stage, stuck,v, w, theta, t_period, c,b,tau,a1,a2,a3):

    t1, t2=t_period
    pi=np.zeros(t2-t1)
    p1=np.zeros(t2-t1)
    p2 = np.full((t2-t1, c), 0.0)
    p3 = np.full((t2-t1, c), 0.0)
    x = np.full((t2-t1, c), 0.0)
    y = np.zeros(t2 - t1)

    for t in range(t1,t2):
        if stage==0 and stuck==False:
            p1[t - t1], p2[t - t1, :], p3[t - t1, :], x[t - t1, :], y[t - t1], pi[t - t1] = subproblem(v, w, theta, t,c, b, tau, a1,a2, a3,till_end=False,time_limit=True, limit=5)
        elif (stage == 1 and stuck == False) or (stage == 0 and stuck == True):
            p1[t - t1], p2[t - t1, :], p3[t - t1, :], x[t - t1, :], y[t - t1], pi[t - t1] = subproblem(v, w, theta, t, c,b,tau,a1,a2,a3, till_end=True, time_limit=True, limit=2*c)
        else:
            p1[t-t1], p2[t-t1,:], p3[t-t1,:], x[t-t1,:], y[t-t1], pi[t - t1] = subproblem(v, w, theta, t,c,b,tau,a1,a2,a3, till_end=True, time_limit=False)

    return [p1, p2, p3, x, y, pi]

def ADP():
    start_time = time.time()
    mp = Model(name='Primal Problem')
    vp = {(t, j): mp.continuous_var(lb=0, name='vp_{}_{}'.format(t, j)) for t in range(T + 1) for j in range(c)}
    thetap = mp.continuous_var_list(T + 1, lb=0, name='thetap')
    wp = mp.continuous_var_list(T + 1, lb=0, name='wp')
    mp.minimize(thetap[T] + mp.sum(vp[T, j] for j in range(c)) + wp[T] * c)

    mp.add_constraint(thetap[0] == 0)
    mp.add_constraint(wp[0] == 0)
    for j in range(c):
        mp.add_constraint(vp[0, j] == 0)
        for t in range(T):
            mp.add_constraint(vp[t, j] <= vp[t + 1, j])
        for t in range(T):
            mp.add_constraint(wp[t] <= wp[t + 1])
            mp.add_constraint(thetap[t] <= thetap[t + 1])

    v = np.full((T + 1, c), 0.0)
    theta = np.zeros(T + 1)
    w = np.zeros(T + 1)
    pi = np.full(T, np.inf)

    Z = 0
    ite = 0
    stage = 0
    time_arr = [0, 0, 0]
    stuck=False
    stuck_time=0

    while sum(pi) > sig * Z or stage <= 1:
        ite = ite + 1
        print('iteration:', ite)
        sub_t_start = time.time()
        if sum(pi) <= sig * Z and stage <= 1:
            time_arr[stage] = time.time() - start_time
            stage += 1

        # Create chunks of the range to distribute across processes
        chunk_size = T // num_processes
        ranges = [(i * chunk_size + 1, (i + 1) * chunk_size + 1) for i in range(num_processes)]
        args = [(stage, stuck,v, w, theta, t_period, c, b, tau, a1, a2, a3) for t_period in ranges]
        # Use a Pool to parallelize the computation
        with Pool(processes=num_processes) as pool:
            results = pool.starmap(execute_subproblem, args)

        p1 = results[0][0]
        p2 = results[0][1]
        p3 = results[0][2]
        x = results[0][3]
        y = results[0][4]
        pi = results[0][5]

        for process in range(1, num_processes):
            p1 = np.concatenate((p1, results[process][0]), axis=0)
            p2 = np.concatenate((p2, results[process][1]), axis=0)
            p3 = np.concatenate((p3, results[process][2]), axis=0)
            x = np.concatenate((x, results[process][3]), axis=0)
            y = np.concatenate((y, results[process][4]), axis=0)
            pi = np.concatenate((pi, results[process][5]), axis=0)

        for t in range(1,T+1):
            print(x[t-1,:], y[t-1])
            print('pi:', pi[t-1])
            if pi[t - 1] <= 0:
                continue
            mp.add_constraint(thetap[t] - thetap[t - 1] + mp.sum((vp[t, j] - vp[t - 1, j]) * x[t-1,j] for j in range(c))
                              + (wp[t] - wp[t - 1]) * (y[t-1] + 1)
                              + p1[t-1] * wp[t - 1]
                              + mp.sum(p2[t-1,j] * (vp[t - 1, j] + wp[t - 1]) for j in range(c))
                              + mp.sum(
                p3[t-1,j] * (vp[t - 1, j] + vp[t - 1, j - (-1) ** (j + 1)] + 2 * wp[t - 1]) for j in range(c))
                              >= p1[t-1] * (r1(p1[t-1], p2[t-1,:], p3[t-1,:]))
                              + sum(p2[t-1,j] * (r2(j, p1[t-1], p2[t-1,:], p3[t-1,:])) for j in range(c))
                              + sum(p3[t-1,j] * (r3(j, p1[t-1], p2[t-1,:], p3[t-1,:])) for j in range(c)))

        sub_t_end = time.time()
        print('sub total time:', sub_t_end - sub_t_start)

        lp_start = time.time()
        mp.solve()
        lp_end = time.time()
        print('lp time:', lp_end - lp_start)
        for t in range(T + 1):
            for j in range(c):
                v[t, j] = vp[t, j].solution_value
            theta = [thetap[t].solution_value for t in range(T + 1)]
            w = [wp[t].solution_value for t in range(T + 1)]
            Z_new = mp.objective_value
        if abs(Z_new-Z)<=0.1:
            stuck_time+=1
        else:
            stuck_time=0
            stuck=False
        Z=Z_new
        if stuck_time>=5:
            stuck=True
        print(theta)
        print(v)
        print(w)
        print(sum(pi))
        print('Revenue:', Z)
    total_time=time.time() - start_time
    print('time first phase:', time_arr[0])
    print('time second phase:', time_arr[1])
    print('time total:', total_time)
    #mp.dump_as_sav(path="C:\\Users\\siqhe\\PycharmProjects\\pythonProject2\\BB", basename='part_3_BB_lp.' + str(c) + 'choice'+str(choice))
    save(v, w, theta, c, choice)

    return Z


sig=0.001
c = 8
choice=1
#print(c)
T = c * 2
num_processes = 8  # Adjust this to the number of cores available
#print("choice",choice)

a1, a2, a3, b, tau = cases.homo_seats(c)

if __name__ == '__main__':
    ADP()
