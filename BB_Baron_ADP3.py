
import time
import numpy as np
import cvxpy as cp
import warnings
# import mosek
import copy
from heapq import *
import itertools
import sys
import ADP_NL_cases as cases

counter = itertools.count()
warnings.simplefilter("error", RuntimeWarning)
eps = 1e-5

class Branch_Bound():
    #for part 3
    def __init__(self, constraints=[], objective=0, x_vars=[], p1_vars=[], p2j_vars=[], p3j_vars=[],
                 y_vars=[], g_vars=[], s_vars=[], bar=-100):
        self.constraints = constraints
        self.objective = objective
        self.x_vars = x_vars
        self.y_vars = y_vars
        self.g_vars=g_vars
        self.s_vars=s_vars
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
                return bestres_r, p10, p20, p30, x0, y0, s0

            if node.is_integral():
                if res>=bestres_r:
                    bestres_r = node.res_r
                    p10 = [v.value for v in prob.var_dict['p1']]
                    p20 = [v.value for v in prob.var_dict['p2j']]
                    p30 = [v.value for v in prob.var_dict['p3j']]
                    x0 = [v.value for v in prob.var_dict['x']]
                    y0 = [v.value for v in prob.var_dict['y']]
                    s0 = [x0[j] * x0[j - (-1) ** (j + 1)] for j in range(len(x0))]
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
                p10 = [v.value for v in node_temp.var_dict['p1']]
                p20 = [v.value for v in node_temp.var_dict['p2j']]
                p30 = [v.value for v in node_temp.var_dict['p3j']]
                x0 = x_vars_r
                y0 = y_vars_r
                s0 = [x0[j]*x0[j-(-1)**(j+1)] for j in range(len(x0))]
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
                return bestres_r, p10, p20, p30, x0, y0, s0
        end=time.time()
        print("Nodes searched: ", nodecount)
        print("time:", end-start)
        print("normal")
        return bestres_r, p10, p20, p30, x0, y0, s0


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
        s0 = [x0[j] * x0[j - (-1) ** (j + 1)] for j in range(len(x0))]
        return res_r, p10, p20, p30, x0, y0, s0


