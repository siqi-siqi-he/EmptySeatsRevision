import docplex
from docplex.mp.model import Model
import Fig7_cases as cases
import math
import time
import numpy as np
import cvxpy as cp
import multiprocessing as multip
import warnings
import copy
from heapq import *
import itertools
import os

counter = itertools.count()

warnings.simplefilter("error", RuntimeWarning)

eps = 1e-5

class BBNode_till_end():
    #part 3
    def __init__(self, vars=[], constraints=[], objective=0, bool_vars=[], p1_vars=[], p2j_vars=[], 
                 func=[], y_vars=[], y_vars_r=[]):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.bool_vars = bool_vars
        self.p1_vars = p1_vars
        self.p2j_vars = p2j_vars
        self.y_vars = y_vars
        self.children = []
        self.bool_vars_r = []
        self.v_vars_r = []
        self.y_vars_r = y_vars_r
        self.res_r = 1e20
        self.func = func

    def buildProblem(self):
        prob = cp.Problem(cp.Maximize(self.objective),
                          self.constraints)
        return prob

    def is_integral(self):
        xall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.bool_vars])
        yall = all([abs(self.y_vars.value - 1) <= 1e-3 or abs(self.y_vars.value - 0) <= 1e-3])
        return xall * yall

    def can_be_integral(self, bool_vars_r, p1, p2j, y):
        value = self.func(p1, p2j, bool_vars_r, y)
        return value

    def branch(self):
        children = []
        for branch in [0, 1]:
            child = copy.deepcopy(self)
            xmin = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[0]
            ymin = (abs(self.y_vars.value - 0.5))
            if xmin < ymin:
                r = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[2]
                child.constraints.append(r == branch)
            else:
                r = self.y_vars
                child.constraints.append(r == branch)
            child.children = []
            children.append(child)
        return children

    def bbsolve1(self):
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
        bestnode = root
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
                end = time.time()
                return bestres_r, bestnode
            bool_vars_r = ([np.round(v.value) for v in node.bool_vars])
            p2j_vars_r = [v.value for v in node.p2j_vars]
            p1_vars_r = np.squeeze([v.value for v in node.p1_vars])
            y_vars_r = np.round(node.y_vars.value)
            res_r = node.can_be_integral(bool_vars_r=bool_vars_r, p1=p1_vars_r, p2j=p2j_vars_r, 
                                         y=y_vars_r)
            if bestres_r > 1e15:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node
            if res<bestres_r:
                continue
            if res_r >= bestres_r:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node

            elif node.is_integral() and res>bestres_r:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node

            new_nodes = node.branch()
            for new_node in new_nodes:
                heappush(heap, (-res, next(counter),
                                new_node))
        end=time.time()
        return bestres_r, bestnode


class BBNode_exist():
    def __init__(self, vars=[], constraints=[], objective=0, bool_vars=[], p1_vars=[], p2j_vars=[], 
                 func=[], y_vars=[], y_vars_r=[]):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.bool_vars = bool_vars
        self.p1_vars = p1_vars
        self.p2j_vars = p2j_vars
        self.y_vars = y_vars
        self.children = []
        self.bool_vars_r = []
        self.v_vars_r = []
        self.y_vars_r = y_vars_r
        self.res_r = 1e20
        self.func = func

    def buildProblem(self):
        prob = cp.Problem(cp.Maximize(self.objective),
                          self.constraints)
        return prob

    def is_integral(self):
        xall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.bool_vars])
        yall = all([abs(self.y_vars.value - 1) <= 1e-3 or abs(self.y_vars.value - 0) <= 1e-3])
        return xall * yall

    def can_be_integral(self, bool_vars_r, p1, p2j, y):
        value = self.func(p1, p2j, bool_vars_r, y)
        return value

    def branch(self):
        children = []
        for b in [0, 1]:
            n1 = copy.deepcopy(self)
            xmin = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[0]
            ymin = (abs(self.y_vars.value - 0.5))
            if xmin < ymin:
                r = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[2]
                n1.constraints.append(r == b)
            else:
                r = self.y_vars
                n1.constraints.append(r == b)
            n1.children = []
            children.append(n1)
        return children

    def bbsolve(self):
        root = self
        try:
            res = root.buildProblem().solve(solver=cp.MOSEK,verbose=False)
        except cp.SolverError as e:
            print(e)
            try:
                res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
            except cp.error.SolverError as e:
                print(e)
                res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)
        heap = [(-res, next(counter), root)]
        bestres_r = 1e20
        bestnode = root
        nodecount = 0
        start=time.time()
        while len(heap) > 0:
            nodecount += 1
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

            if prob.status in ["infeasible","unbounded"]:
                end=time.time()
                return bestres_r, bestnode
            bool_vars_r = ([np.round(v.value) for v in node.bool_vars])
            p2j_vars_r = [v.value for v in node.p2j_vars]
            p1_vars_r = np.squeeze([v.value for v in node.p1_vars])
            y_vars_r = np.round(node.y_vars.value)
            res_r = node.can_be_integral(bool_vars_r=bool_vars_r, p1=p1_vars_r, p2j=p2j_vars_r, 
                                         y=y_vars_r)
            if bestres_r > 1e15:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node
                end = time.time()
                return bestres_r, bestnode

            if res<bestres_r:
                continue
            if res_r >= bestres_r:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node

            elif node.is_integral():
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                node.y_vars_r = y_vars_r
                bestnode = node

            new_nodes = node.branch()
            for new_node in new_nodes:
                heappush(heap, (-res, next(counter),
                                new_node))
        end=time.time()
        return bestres_r, bestnode




def kl_div(x, y):
    return x * (math.log(x) - math.log(y)) - x + y



def r1(p1, p2):
    if p1 <= 0:
        result = 0
        print('r1 is bad', p1)
    else:
        try:
            result = (a1 -a0/tau[0]) / b[0] + 1 / b[0] * 1 / tau[0] * (math.log(1 - p1 - sum(p2)) - math.log(p1))
        except ValueError as e:
            # Print the variable causing the ValueError
            print(f"ValueError: {e}")
            print("Printing variable values:")
            print(f"p1: {p1}")
            print(f"p2: {p2}")
            print(f"tau[2]: {tau[2]}")
            result = e
    return result

def r2(j, p1, p2):
    try:
        if p2[j] <= 0:
            print('r2 is bad', p2[j])
            result = 0
        elif sum(p2) <= 0:
            for i in range(len(p2)):
                if p2[i] <= 0:
                    p2[i] = 0
            result = (a2[j] -a0/tau[1])/ b[1] + 1 / b[1] * (math.log(1 - p1 - sum(p2) ) - math.log(p2[j])) + 1 / b[1] * (
                    1 - tau[1]) / tau[1] * (math.log(1 - p1 - sum(p2) ) - math.log(sum(p2)))
        else:
            result = (a2[j] -a0/tau[1]) / b[1] + 1 / b[1] * (math.log(1 - p1 - sum(p2) ) - math.log(p2[j])) + 1 / b[1] * (
                        1 - tau[1]) / tau[1] * (math.log(1 - p1 - sum(p2) ) - math.log(sum(p2)))
    except ValueError as e:
        # Print the variable causing the ValueError
        print(f"ValueError: {e}")
        print("Printing variable values:")
        print(f"a2[j]: {a2[j]}")
        print(f"b[2]: {b[2]}")
        print(f"p1: {p1}")
        print(f"p2: {p2}")
        print(f"tau[1]: {tau[1]}")
        print(f"j: {j}")
    return result

def obj_func(v, w, theta, t, p10, p20,, x, y):
    obj_value = p10 * r1(p10, p20) - w[t - 1] * p10 + sum(
        [p20[j] * (r2(j, p10, p20) - v[t - 1, j] - w[t - 1]) for j in range(c)]) \
                - theta[t] + theta[t - 1] - (w[t] - w[t - 1]) * (y + 1) - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
    return obj_value

def subproblem(v, w, theta, t, type):
    # declare/define the decision variables of the subproblem
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p2j = cp.Variable(c, name="p2j")
    x = cp.Variable(c, name="x")
    y = cp.Variable(1, name="y")
    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + (a1-a0/tau[0]) / b[0] * p1 - w[t - 1] * p1 \
                   + cp.sum(
        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + (a2[j]-a0/tau[1]) / b[1] * p2j[j] - (v[t - 1, j] + w[t - 1]) * p2j[j]
         for j in range(c)]) \
                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                   - theta[t] + theta[t - 1] - (w[t] - w[t - 1]) * (y + 1) - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
    constraints = [p1 <= 1,
                   p1 >= eps,
                   p0 >= eps,
                   p2 >= eps,
                   p2 <= 1,
                   p1 + p0 + p2 == 1,
                   p2 == cp.sum(p2j),
                   p0 >= 0,
                   y <= 1,
                   y >= 0,
                   x >= 0,
                   x <= 1,
                   p1 <= y + 1,
                   p2j <= y + 1,
                   p2j <= 1,
                   p2j >= eps / c,
                   p2j <= x + eps / c,
                   cp.sum(x)>=y+1]
    bool_vars = [x[i] for i in range(c)]
    y_vars = {y}
    vars = {p0[0], p2[0]}
    func = lambda p10, p20, p30, x, y: obj_func(v, w, theta, t, p10, p20, p30, x, y)
    if type==0:
        root = BBNode_exist(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars, p1_vars={p1},
                          p2j_vars={p2j[j] for j in range(c)}, func=func, y_vars=y)
        res, sol = root.bbsolve()
    if type==1:
        root = BBNode_till_end(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars, p1_vars={p1},
                          p2j_vars={p2j[j] for j in range(c)}, func=func, y_vars=y)
        res, sol = root.bbsolve1()
    p10 = np.squeeze([v.value for v in sol.p1_vars])
    p20 = [v.value for v in sol.p2j_vars]
    if res <= 0:
        x0 = [v.value for v in sol.bool_vars]
        y0 = sol.y_vars[0].value
    else:
        x0 = sol.bool_vars_r
        y0 = sol.y_vars_r[0]
    obj_value = res
    return p10, p20, x0, y0, obj_value

def save(v, w, theta):

    folder_path = "wo3_ADP"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "v_value_choice" + str(preference) + "a0" + str(a0) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, v)
    file_name = "w_value_choice" + str(preference) + "a0" + str(a0) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, w)
    file_name = "theta_value_choice" + str(preference) + "a0" + str(a0) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, theta)

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
    while sum(pi) > sig * Z:
        sub_t_start = time.time()
        for t in range(1, T + 1):
            p1, p2, x, y, pi[t - 1] = subproblem(v, w, theta, t, 0)

            r10=r1(p1, p2)
            r20=[r2(j, p1, p2) for j in range(c)]
            flag = 0
            for i in range(c):
                if r20[i] < 0:
                    flag = 1
            if r10 < 0 or flag == 1:
                raise Exception("negative prices!", r10, r20, r30)
            if pi[t - 1] <= 0:
                continue
            mp.add_constraint(thetap[t] - thetap[t - 1] + mp.sum((vp[t, j] - vp[t - 1, j]) * x[j] for j in range(c))
                              + (wp[t] - wp[t - 1]) * (y + 1)
                              + p1 * wp[t - 1]
                              + mp.sum(p2[j] * (vp[t - 1, j] + wp[t - 1]) for j in range(c))
                              >= p1 * (r1(p1, p2))
                              + sum(p2[j] * (r2(j, p1, p2)) for j in range(c)))
        sub_t_end = time.time()

        lp_start = time.time()
        mp.solve()
        lp_end = time.time()


        for t in range(T + 1):
            for j in range(c):
                v[t, j] = vp[t, j].solution_value
            theta = [thetap[t].solution_value for t in range(T + 1)]
            w = [wp[t].solution_value for t in range(T + 1)]
            Z = mp.objective_value

    end_time = time.time()
    elapsed_time = end_time - start_time
    pi = np.full(T, np.inf)

    time_t_s=time.time()
    while sum(pi) > sig * Z:
        sub_start = time.time()
        for t in range(1, T + 1):
            sub_t_start = time.time()
            p1, p2,x, y, pi[t - 1] = subproblem(v, w, theta, t,1)
            r10 = r1(p1, p2)
            r20 = [r2(j, p1, p2) for j in range(c)]
            flag=0
            for i in range(c):
                if r20[i]<0 :
                    flag=1
            if r10<0 or flag==1:
                raise Exception("negative prices!",r10)
            if pi[t - 1] <= 0:
                continue
            mp.add_constraint(thetap[t] - thetap[t - 1] + mp.sum((vp[t, j] - vp[t - 1, j]) * x[j] for j in range(c))
                              + (wp[t] - wp[t - 1]) * (y + 1)
                              + p1 * wp[t - 1]
                              + mp.sum(p2[j] * (vp[t - 1, j] + wp[t - 1]) for j in range(c))
                              >= p1 * (r1(p1, p2))
                              + sum(p2[j] * (r2(j, p1, p2)) for j in range(c)))
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
    time_t_e=time.time()
    #print('total part 3', time_t_e-time_t_s)
    save(v, w, theta)
    return Z

sig=0.001
eps = 1e-5
results=np.zeros((3,21))
c=8
T=16
a1, a2, a3, b, tau = cases.homo_seats(c)
folder_path="AveragePrices"
for preference in range(1,2):
    #high, medium, low
    #+-0.2
    for strength in range(21):
        a0=strength*0.5-5
        if preference==0:
            file_name="a1_a0_"+str(a0)+"_a3decr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a1=np.loadtxt(file_path)
            file_name="a2_a0_"+str(a0)+"_a3decr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a2=np.loadtxt(file_path)
            a3=[0.4]*8
        if preference==1:
            a1=0.2
            a3=[0.6]*8
            a2=[0.4]*8
        if preference==2:
            file_name="a1_a0_"+str(a0)+"_a3incr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a1=np.loadtxt(file_path)
            file_name="a2_a0_"+str(a0)+"_a3incr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a2=np.loadtxt(file_path)
            a3=[0.8]*8
        print(a0)
        results[preference][strength]=ADP()
        print(results[preference][strength])
    print(results)