import cvxpy as cp


x = cp.Variable(1, name="x")

objective_cp1= x*x
objective_cp2= cp.multiply(x,x)
objective_cp3= cp.square(x)

print(objective_cp1.curvature)
print(objective_cp2.curvature)
print(objective_cp3.curvature)