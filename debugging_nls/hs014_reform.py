# it was useful to have these evaluations laying around. maybe it will be again
f = lambda x: (x[0]-2)**2 + (x[1]-1)**2
c1 = lambda x: .25*(x[0])**2 + x[1]**2 - 1
c2 = lambda x: x[0] - 2*x[1] + 1

df = lambda x: [2*(x[0]-2), 2*(x[1] - 1)]
dc1 = lambda x: [2*.25*x[0], 2*x[1]]
dc2 = lambda x: [1,-2]

x0 = [2,2]
print(f(x0))
print(c1(x0))
print(c2(x0))

print(df(x0))
print(dc1(x0))
print(dc2(x0))

# print

from math import sqrt
norm = lambda le,li: sqrt( (le + li)**2 + (4*li-2*le)**2)
print(norm(0.0, 0.0))
