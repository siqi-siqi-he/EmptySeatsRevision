import docplex
from docplex.mp.model import Model
import new_cases as cases
import math
import time
import numpy as np
import cvxpy as cp
import warnings
# import mosek
import copy
from heapq import *
import itertools
import sys
import os
import pandas as pd
from multiprocessing import Pool, Manager

counter = itertools.count()
#warnings.simplefilter("error", RuntimeWarning)
eps = 1e-5

#This is the branch and bound for optimal solution
class BBTreeNode1():
    #for part 3
    def __init__(self, vars=[], constraints=[], objective=0, bool_vars=[], p1_vars=[], p2j_vars=[], p3j_vars=[],
                 func=[], y_vars=[], y_vars_r=[], capacity=[]):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.bool_vars = bool_vars
        self.p1_vars = p1_vars
        self.p2j_vars = p2j_vars
        self.p3j_vars = p3j_vars
        self.y_vars = y_vars
        self.children = []
        self.bool_vars_r = []
        self.v_vars_r = []
        self.y_vars_r = y_vars_r
        self.res_r = 1e20
        self.func = func
        self.c=capacity

    def buildProblem(self):
        prob = cp.Problem(cp.Maximize(self.objective),
                          self.constraints)  # i put Minimize, just so you know that I'm assuming it
        return prob

    def is_integral(self):
        xall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.bool_vars])
        yall = all([abs(self.y_vars.value - 1) <= 1e-3 or abs(self.y_vars.value - 0) <= 1e-3])
        return xall * yall

    def can_be_integral(self, bool_vars_r, p1, p2j, p3j, y):
        value = self.func(p1, p2j, p3j, bool_vars_r, y)
        return value

    def branch(self):
        children = []
        for branch in [0, 1]:
            child = copy.deepcopy(self)  # yeesh. Not good performance wise, but is simple implementation-wise
            xmin = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[0]
            ymin = (abs(self.y_vars.value - 0.5))
            if xmin < ymin:
                r = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[2]
                child.constraints.append(r == b)
            #                n1.bool_vars.remove(r)  # remove binary constraint from bool var set
            else:
                r = self.y_vars
                child.constraints.append(r == branch)  # add in the new binary constraint
            #                n1.y_vars.remove(r)  # remove binary constraint from bool var set

            child.children = []
            #            n1.vars.add(v)  # and add it into var set for later inspection of answer
            children.append(child)
        return children

    def bbsolve1(self, ite):
        root = self
        #ite_max=78
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
        bestnode = root  # initialize bestnode to the root
        #print(heap)
        nodecount = 0
        start=time.time()
        while len(heap) > 0:
            nodecount += 1  # for statistics
            # print("Heap Size: ", len(heap))
            _, _, node = heappop(heap)
            prob = node.buildProblem()
            #            print([v.value for v in node.bool_vars])

            try:
                res = prob.solve(solver=cp.MOSEK,verbose=False)
            except cp.error.SolverError as e:
                print(e)
                try:
                    res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
                except cp.error.SolverError as e:
                    print(e)
                    res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)

            if prob.status == "infeasible":
                end = time.time()

                print("Nodes searched: ", nodecount)
                print("time:", end - start)
                return bestres_r, bestnode
            #            print([v.value for v in node.bool_vars])
            bool_vars_r = ([np.round(v.value) for v in node.bool_vars])

            p2j_vars_r = [v.value for v in node.p2j_vars]
            p3j_vars_r = [v.value for v in node.p3j_vars]
            p1_vars_r = np.squeeze([v.value for v in node.p1_vars])
            y_vars_r = np.round(node.y_vars.value)

            res_r = node.can_be_integral(bool_vars_r=bool_vars_r, p1=p1_vars_r, p2j=p2j_vars_r, p3j=p3j_vars_r,
                                         y=y_vars_r)


            if bestres_r > 1e15:
                bestres_r = res_r
            if res_r > 0 and res_r >= bestres_r:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node

            # print("Result: ", res)
            if prob.status not in ["infeasible", "unbounded"]:
                if res < 0 + 1e-5:  # even the relaxed problem is feasible. forget about this branch then
                    # print("Relaxed Problem Feasible. Killing this branch.")
                    pass
                if res < bestres_r - 1e-3:  # even the relaxed problem sucks. forget about this branch then
                    # print("Relaxed Problem Stinks. Killing this branch.")
                    pass
                elif node.is_integral():

                    bestres_r = res_r
                    node.res_r = res_r
                    node.bool_vars_r = bool_vars_r
                    node.y_vars_r = y_vars_r
                    bestnode = node

                    end = time.time()

                    print("Nodes searched: ", nodecount)
                    print("time:", end - start)
                    return bestres_r, bestnode
                    if res > bestres_r: # if a valid solution then this is the new best
                        # print("New Best Integral solution.")
                        bestres_r = res
                        node.res_r = res
                        node.bool_vars_r = [node.bool_vars[i].value for i in range(self.c)]
                        node.y_vars_r = node.y_vars.value
                        bestnode = node

                else:  # otherwise, we're unsure if this branch holds promise. Maybe it can't actually achieve this lower bound. So branch into it
                    new_nodes = node.branch()
                    for new_node in new_nodes:
                        heappush(heap, (-res, next(counter),
                                        new_node))  # using counter to avoid possible comparisons between nodes. It tie breaks
            else:
                print("Problem status:", prob.status)
        end=time.time()

        print("Nodes searched: ", nodecount)
        print("time:", end-start)
        return bestres_r, bestnode

#This is the branch and bound for existing solution
class BBTreeNode():
    def __init__(self, vars=[], constraints=[], objective=0, bool_vars=[], p1_vars=[], p2j_vars=[], p3j_vars=[],
                 func=[], y_vars=[], y_vars_r=[],capacity=[]):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.bool_vars = bool_vars
        self.p1_vars = p1_vars
        self.p2j_vars = p2j_vars
        self.p3j_vars = p3j_vars
        self.y_vars = y_vars
        self.children = []
        self.bool_vars_r = []
        self.v_vars_r = []
        self.y_vars_r = y_vars_r
        self.res_r = 1e20
        self.func = func
        self.c=capacity

    def buildProblem(self):
        prob = cp.Problem(cp.Maximize(self.objective),
                          self.constraints)
        return prob

    def is_integral(self):
        #check whether the solution is integral
        xall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.bool_vars])
        yall = all([abs(self.y_vars.value - 1) <= 1e-3 or abs(self.y_vars.value - 0) <= 1e-3])
        return xall * yall

    def can_be_integral(self, bool_vars_r, p1, p2j, p3j, y):
        value = self.func(p1, p2j, p3j, bool_vars_r, y)
        return value

    def branch(self):
        children = []
        for branch in [0, 1]:
            child = copy.deepcopy(self)
            xmin = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[0]
            ymin = (abs(self.y_vars.value - 0.5))
            if xmin < ymin:
                r = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[2]
                child.constraints.append(r == b)
            else:
                r = self.y_vars
                child.constraints.append(r == branch)  # add in the new binary constraint
            child.children = []
            children.append(child)
        return children

    def bbsolve(self, ite):
        root = self
        ite_max=3000

        try:
            #solve the root problem
            res = root.buildProblem().solve(solver=cp.MOSEK,verbose=False)
        except cp.SolverError as e:
            print(e)
            try:
                res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
            except cp.error.SolverError as e:
                print(e)
                res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)

        #define a priority queue
        heap = [(-res, next(counter), root)]
        bestres_r = 1e20
        bestnode = root  # initialize bestnode to the root
        #print(heap)
        nodecount = 0
        start=time.time()
        while len(heap) > 0:
            nodecount += 1
            # print("Heap Size: ", len(heap))
            #pop a node from the priority queue, then solve it
            _, _, node = heappop(heap)
            prob = node.buildProblem()

            try:
                res = prob.solve(solver=cp.MOSEK,verbose=False)
            except cp.SolverError as e:
                print(e)
                try:
                    res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
                except cp.error.SolverError as e:
                    print(e)
                    res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)

            if prob.status == "infeasible":
                return bestres_r, bestnode

            bool_vars_r = ([np.round(v.value) for v in node.bool_vars])
            p2j_vars_r = [v.value for v in node.p2j_vars]
            p3j_vars_r = [v.value for v in node.p3j_vars]
            p1_vars_r = np.squeeze([v.value for v in node.p1_vars])
            y_vars_r = np.round(node.y_vars.value)
            res_r = node.can_be_integral(bool_vars_r=bool_vars_r, p1=p1_vars_r, p2j=p2j_vars_r, p3j=p3j_vars_r,
                                         y=y_vars_r)
            if ite<ite_max:
                if res_r>0:
                    bestres_r = res_r
                    node.res_r = res_r
                    node.bool_vars_r = bool_vars_r
                    node.y_vars_r = y_vars_r
                    bestnode = node

                    return bestres_r, bestnode
            if bestres_r > 1e15:
                bestres_r = res_r
            if res_r > 0 and res_r >= bestres_r:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node

            # print("Result: ", res)
            if prob.status not in ["infeasible", "unbounded"]:
                if res < 0 + 1e-5:  # even the relaxed problem is feasible. forget about this branch then
                    # print("Relaxed Problem Feasible. Killing this branch.")
                    pass
                if res < bestres_r - 1e-3:  # even the relaxed problem sucks. forget about this branch then
                    # print("Relaxed Problem Stinks. Killing this branch.")
                    pass
                elif node.is_integral():
                    if ite<ite_max:
                        bestres_r = res_r
                        node.res_r = res_r
                        node.bool_vars_r = bool_vars_r
                        node.y_vars_r = y_vars_r
                        bestnode = node

                        return bestres_r, bestnode
                    if res > bestres_r: # if a valid solution then this is the new best
                        # print("New Best Integral solution.")
                        bestres_r = res
                        node.res_r = res
                        node.bool_vars_r = [node.bool_vars[i].value for i in range(self.c)]
                        node.y_vars_r = node.y_vars.value
                        bestnode = node

                else:  # otherwise, we're unsure if this branch holds promise. Maybe it can't actually achieve this lower bound. So branch into it
                    new_nodes = node.branch()
                    for new_node in new_nodes:
                        heappush(heap, (-res, next(counter),
                                        new_node))  # using counter to avoid possible comparisons between nodes. It tie breaks
            else:
                print("Problem status:", prob.status)
            end=time.time()
            if end-start>=60:
                break
        print("Nodes searched: ", nodecount)

        return bestres_r, bestnode

def UP_t(c, a1, a2, a3, b, tau):
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def UP_1(c, a1, a2, a3, b, tau):
    result = math.exp(tau[0] * (a1)) / (math.exp(tau[0] * (a1)) + 1)
    return result

def kl_div(x, y):
    # this function corresponds to kl_div(x,y) in the function library provided by cvxpy
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
    #the price for product 1 given the purchase probability p1, p2, p3
    if p1 <= 0:
        result = 0
        print('r1 is bad', p1)
    else:
        result = a1 / b[0] + 1 / b[0] * 1 / tau[0] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p1))
    return result

def r2(j, p1, p2, p3):
    # the price array for product 2 given the purchase probability p1, p2, p3
    if p2[j] <= 0:
        print('r2 is bad', p2[j])
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
    # the price array for product 3 given the purchase probability p1, p2, p3
    if p3[j] <= 0:
        result = 0
        print('r3 is bad', p3[j])
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

def obj_func(v, w, theta, t, p10, p20, p30, x, y,c):
    #objective function of the subproblem, in case you need it for debugging, etc
    obj_value = p10 * r1(p10, p20, p30) - w[t - 1] * p10 + sum(
        [p20[j] * (r2(j, p10, p20, p30) - v[t - 1, j] - w[t - 1]) for j in range(c)]) \
                + sum(
        [p30[j] * (r3(j, p10, p20, p30) - v[t - 1, j] - v[t - 1, j - (-1) ** (j + 1)] - 2 * w[t - 1]) for j in
         range(c)]) \
                - theta[t] + theta[t - 1] - (w[t] - w[t - 1]) * (y + 1) - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
    return obj_value

def subproblem(v, w, theta, t, ite,type,c):
    # the subproblem function finds the violating constraints in ADP
    # if type=0, it finds the existing violating constraints
    # if type=1, it finds the most violating constraints
    # declare/define the decision variables of the subproblem
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    x = cp.Variable(c, name="x")
    y = cp.Variable(1, name="y")

    #define objective by disciplined convex programming
    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 - w[t - 1] * p1 \
                   + cp.sum(
        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j])
         + a2[j] / b[1] * p2j[j]
         - (v[t - 1, j] + w[t - 1]) * p2j[j]
         for j in range(c)]) \
                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                   + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (
                v[t - 1, j] + v[t - 1, j - (-1) ** (j + 1)] + 2 * w[t - 1]) * p3j[j] for j in range(c)]) \
                   + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) \
                   - theta[t] + theta[t - 1] - (w[t] - w[t - 1]) * (y + 1) - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))

    #Constructing constraints
    constraints = [p1 <= 1,
                   p1 >= eps,
                   p0 >= eps,
                   p2 >= eps,
                   p2 <= 1,
                   p3 >= eps,
                   p3 <= 1,
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
                   p3j <= y + eps / c,  # p3j <= math.floor(y / 2) + eps / c,
                   p2j <= 1,
                   p3j <= 1,
                   p2j >= eps / c,
                   p3j >= eps / c,
                   p2j <= x + eps / c,
                   p3j <= x + eps / c,
                   cp.sum(x)>=y+1]

    for j in range(c):
        constraints = constraints + [p3j[j] <= x[j - (-1) ** (j + 1)] + eps / c]

    bool_vars = [x[i] for i in range(c)]
    y_vars = {y}
    vars = {p0[0], p2[0], p3[0]}
    func = lambda p10, p20, p30, x, y: obj_func(v, w, theta, t, p10, p20, p30, x, y,c)
    if type==0:
        #bnb for existing solution
        root = BBTreeNode(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars, p1_vars={p1},
                          p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)}, func=func, y_vars=y, capacity=c)
        res, sol = root.bbsolve(ite)
    if type==1:
        #bnb for optimal solution
        root = BBTreeNode1(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars, p1_vars={p1},
                          p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)}, func=func, y_vars=y, capacity=c)
        res, sol = root.bbsolve1(ite)

    p10 = np.squeeze([v.value for v in sol.p1_vars])
    p20 = [v.value for v in sol.p2j_vars]
    p30 = [v.value for v in sol.p3j_vars]
    if res <= 0:
        x0 = [v.value for v in sol.bool_vars]
        y0 = sol.y_vars[0].value

    else:
        x0 = sol.bool_vars_r
        y0 = sol.y_vars_r[0]
    obj_value = res

    return p10, p20, p30, x0, y0, obj_value

def execute_subproblem(v, w, theta, t_period, ite,type,c):
    #this function is invented only for parallel computing
    t1, t2=t_period
    pi=np.zeros(t2-t1)
    p1=np.zeros(t2-t1)
    p2 = np.full((t2-t1, c), 0.0)
    p3 = np.full((t2-t1, c), 0.0)
    x = np.full((t2-t1, c), 0.0)
    y = np.zeros(t2 - t1)

    for t in range(t1,t2):
        p1[t - t1], p2[t - t1, :], p3[t - t1, :], x[t - t1, :], y[t - t1], pi[t - t1] = subproblem(v, w, theta, t, ite, type,c)

    return [p1, p2, p3, x, y, pi]

def save(v, w, theta, c, choice):
    #the results of affine function are saved
    folder_path = "ADP"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, v)
    file_name = "w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, w)
    file_name = "theta_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, theta)

def ADP(c):

    start_time = time.time()
    mp = Model(name='Primal Problem')
    #you can change the method when solving the lp problem
    #mp.parameters.lpmethod=2

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

    while sum(pi) > sig * Z:
        ite = ite + 1
        print('iteration:', ite)
        sub_t_start = time.time()

        chunk_size = T // num_processes
        ranges = [(i * chunk_size + 1, (i + 1) * chunk_size + 1) for i in range(num_processes)]
        args = [(v, w, theta, t_period, ite, 0,c) for t_period in ranges]
        #print("before para")
        with Pool(processes=num_processes) as pool:
            results = pool.starmap(execute_subproblem, args)
        #print("after para")
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

        for t in range(1, T + 1):
            print(x[t - 1, :], y[t - 1])
            print('pi:', pi[t - 1])
            if pi[t - 1] <= 0:
                continue

            mp.add_constraint(
                thetap[t] - thetap[t - 1] + mp.sum((vp[t, j] - vp[t - 1, j]) * x[t - 1, j] for j in range(c))
                + (wp[t] - wp[t - 1]) * (y[t - 1] + 1)
                + p1[t - 1] * wp[t - 1]
                + mp.sum(p2[t - 1, j] * (vp[t - 1, j] + wp[t - 1]) for j in range(c))
                + mp.sum(
                    p3[t - 1, j] * (vp[t - 1, j] + vp[t - 1, j - (-1) ** (j + 1)] + 2 * wp[t - 1]) for j in range(c))
                >= p1[t - 1] * (r1(p1[t - 1], p2[t - 1, :], p3[t - 1, :]))
                + sum(p2[t - 1, j] * (r2(j, p1[t - 1], p2[t - 1, :], p3[t - 1, :])) for j in range(c))
                + sum(p3[t - 1, j] * (r3(j, p1[t - 1], p2[t - 1, :], p3[t - 1, :])) for j in range(c)))

        sub_t_end = time.time()
        #print('sub total time:', sub_t_end - sub_t_start)

        lp_start = time.time()
        mp.solve()
        lp_end = time.time()
        #print('lp time:', lp_end - lp_start)

        for t in range(T + 1):
            for j in range(c):
                v[t, j] = vp[t, j].solution_value
            theta = [thetap[t].solution_value for t in range(T + 1)]
            w = [wp[t].solution_value for t in range(T + 1)]
            Z = mp.objective_value
        print(theta)
        print(v)
        print(w)
        print(sum(pi))
        print('Revenue:', Z)
    end_time = time.time()
    elapsed_time = end_time - start_time
    #print('time in total:', elapsed_time)

    ite_before_3=ite

    #the second round
    pi = np.full(T, np.inf)
    #print('optimality solving starts')
    time_t_s=time.time()
    while sum(pi) > sig * Z:
        ite = ite + 1
        print('iteration:', ite)
        sub_start = time.time()
        chunk_size = T // num_processes
        ranges = [(i * chunk_size + 1, (i + 1) * chunk_size + 1) for i in range(num_processes)]
        args = [(v, w, theta, t_period, ite, 1,c) for t_period in ranges]
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

        for t in range(1, T + 1):
            sub_t_start = time.time()
            print(x[t - 1, :], y[t - 1])
            print('pi:', pi[t - 1])
            if pi[t - 1] <= 0:
                continue

            mp.add_constraint(
                thetap[t] - thetap[t - 1] + mp.sum((vp[t, j] - vp[t - 1, j]) * x[t - 1, j] for j in range(c))
                + (wp[t] - wp[t - 1]) * (y[t - 1] + 1)
                + p1[t - 1] * wp[t - 1]
                + mp.sum(p2[t - 1, j] * (vp[t - 1, j] + wp[t - 1]) for j in range(c))
                + mp.sum(
                    p3[t - 1, j] * (vp[t - 1, j] + vp[t - 1, j - (-1) ** (j + 1)] + 2 * wp[t - 1]) for j in range(c))
                >= p1[t - 1] * (r1(p1[t - 1], p2[t - 1, :], p3[t - 1, :]))
                + sum(p2[t - 1, j] * (r2(j, p1[t - 1], p2[t - 1, :], p3[t - 1, :])) for j in range(c))
                + sum(p3[t - 1, j] * (r3(j, p1[t - 1], p2[t - 1, :], p3[t - 1, :])) for j in range(c)))

            sub_t_end = time.time()
            #print('sub single time:', sub_t_end - sub_t_start)
        sub_end=time.time()
        #print('sub total time', sub_end-sub_start)
        lp_start = time.time()
        mp.solve()
        lp_end = time.time()
        #print('lp time:', lp_end - lp_start)
        for t in range(T + 1):
            for j in range(c):
                v[t, j] = vp[t, j].solution_value
            theta = [thetap[t].solution_value for t in range(T + 1)]
            w = [wp[t].solution_value for t in range(T + 1)]
            Z = mp.objective_value
        if __name__ == "__main__":
            print(theta)
            print(v)
            print(w)
            print(sum(pi))
            print('Revenue:', Z)
    time_t_e=time.time()
    total_3=time_t_e-time_t_s
    #print('total part for solving to optimality', total_3)
    save(v, w, theta, c, choice)
    print('capacity', c, 'choice', choice, 'Revenue', Z, 'time 0', elapsed_time, 'time 1', total_3, 'ite_before', ite_before_3, 'ite_3', ite,)


    return Z, time_t_e-start_time



sig=0.001
#set the number of processors for parallel computing
num_processes=16
#set the capacity c
c=48
#set the selling horizon days
T = c * 2
#set choice
choice=1
if choice==1:
    a1, a2, a3, b, tau = cases.homo_seats(c)
elif choice==2:
    a1, a2, a3, b, tau = cases.aw_seats(c)
elif choice==3:
    a1, a2, a3, b, tau = cases.incre_seats(c)

folder_path = "results"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
file_name = "runtime_ADP" + str(choice) + "capacity.csv"
file_path_time = f"{folder_path}/{file_name}"
try:
    df_time = pd.read_csv(file_path_time, index_col=0)
except FileNotFoundError:
    df_time = pd.DataFrame(columns=["i", "Result"]).set_index("i")

file_name = "UB_ADP" + str(choice) + "capacity.csv"
file_path_UB = f"{folder_path}/{file_name}"
try:
    df_UB = pd.read_csv(file_path_UB, index_col=0)
except FileNotFoundError:
    df_UB = pd.DataFrame(columns=["i", "Result"]).set_index("i")

if __name__ == "__main__":
    UB, Time=ADP(c)
    df_time.loc[int(c/8)] = Time
    df_UB.loc[int(c/8)]=UB
    df_time.to_csv(file_path_time)
    df_UB.to_csv(file_path_UB)
