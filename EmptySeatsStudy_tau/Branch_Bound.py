import copy
from heapq import *
import itertools
import sys
import time
import numpy as np
import cvxpy as cp

counter = itertools.count()


class BBNode():
    #search till the end
    def __init__(self, vars=[], constraints=[], objective=0, bool_vars=[], p0_vars=[],p1_vars=[], p2j_vars=[], p3j_vars=[],
                 func=[]):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.bool_vars = bool_vars
        self.p0_vars = p0_vars
        self.p1_vars = p1_vars
        self.p2j_vars = p2j_vars
        self.p3j_vars = p3j_vars
        self.children = []
        self.bool_vars_r = []
        self.v_vars_r = []
        self.res_r = 1e20
        self.func = func

    def buildProblem(self):
        prob = cp.Problem(cp.Maximize(self.objective),
                          self.constraints)
        return prob

    def is_integral(self):
        xall = all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.bool_vars])
        return xall

    def can_be_integral(self, bool_vars_r, p0,p1, p2j, p3j):
        value = self.func(p0,p1, p2j, p3j, bool_vars_r)
        return value

    def branch(self):
        children = []
        for branch in [0, 1]:
            child = copy.deepcopy(self)
            r = min([(abs(v.value - 0.5), i, v) for i, v in enumerate(self.bool_vars)])[2]
            child.constraints.append(r == branch)
            child.children = []
            children.append(child)
        return children

    def bbsolve(self):
        root = self
        try:
            res = root.buildProblem().solve(solver=cp.MOSEK, verbose=False)
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
                res = prob.solve(solver=cp.MOSEK, verbose=False)
            except cp.error.SolverError as e:
                print(e)
                try:
                    res = root.buildProblem().solve(solver=cp.SCS, verbose=False)
                except cp.error.SolverError as e:
                    print(e)
                    res = root.buildProblem().solve(solver=cp.ECOS, verbose=False)

            if prob.status in ["infeasible", "unbounded"]:
                end = time.time()
                #print("Nodes searched: ", nodecount)
                #print("time:", end - start)
                return bestres_r, bestnode
            bool_vars_r = ([np.round(v.value) for v in node.bool_vars])
            p2j_vars_r = [v.value for v in node.p2j_vars]
            p3j_vars_r = [v.value for v in node.p3j_vars]
            p1_vars_r = np.squeeze([v.value for v in node.p1_vars])
            p0_vars_r = np.squeeze([v.value for v in node.p0_vars])
            res_r = node.can_be_integral(bool_vars_r=bool_vars_r, p0=p0_vars_r,p1=p1_vars_r, p2j=p2j_vars_r, p3j=p3j_vars_r)

            if bestres_r > 1e15:
                bestres_r = res_r
            if res<bestres_r:
                continue
            if res_r >= bestres_r:
                bestres_r = res_r
                node.res_r = res_r
                node.bool_vars_r = bool_vars_r
                bestnode = node

            elif node.is_integral() and res>bestres_r:
                bestres_r = res
                node.res_r = res
                node.bool_vars_r = bool_vars_r
                bestnode = node

            new_nodes = node.branch()
            for new_node in new_nodes:
                heappush(heap, (-res, next(counter),
                                new_node))
        end=time.time()
        #print("Nodes searched: ", nodecount)
        #print("time:", end-start)
        return bestres_r, bestnode
