from fenics import *
from matplotlib import pyplot as plt
import numpy as np
from mshr import *

from Droplet_Numerics.meshes_initialData_postProcessing.Set_2.model_parameters    import *

mesh    = Mesh('/content/Droplet_Numerics/meshes_initialData_postProcessing/mesh.xml')  

P2        = VectorElement("P", mesh.ufl_cell(), 1)
P1        = FiniteElement("P", mesh.ufl_cell(), 1)
R         = FiniteElement("P", mesh.ufl_cell(), 1)    

TH        = MixedElement([P2,R]) 
W         = FunctionSpace(mesh, TH)
S         = FunctionSpace(mesh, P1) 

psi11 = Expression(("tanh(ff*(h0-x[1])/eps)"),eps=eps,h0=h0,ff=interface_factor,degree=2)
psi22 = Expression(("-1 + 2*(1+tanh(ff*(x[1]/h0-h0)/eps))*(1+tanh(ff/eps*(h1-pow(pow((x[0])/h0-H1/2,2.0)+pow(x[1]/h0-h0,2.0),0.5))))/4"),eps=eps,h0=h0,h1=h1,H1=H1,ff=interface_factor,degree=2)    
    

initial = Expression(("0","0","0"), degree = 2) 
old_q   = interpolate(initial, W)
 
bc, W, S, P2, mesh, psi1, psi2  = refine_mesh_pureDef(mesh,W,S,old_q,eps,h0,h1,interface_factor,H1)
old_q                           = interpolate(initial, W)
bc, W, S, P2, mesh, psi1, psi2  = refine_mesh_pureDef(mesh,W,S,old_q,eps,h0,h1,interface_factor,H1)
old_q                           = interpolate(initial, W)
bc, W, S, P2, mesh, psi1, psi2  = refine_mesh_pureDef(mesh,W,S,old_q,eps,h0,h1,interface_factor,H1)
old_q                           = interpolate(initial, W)
bc, W, S, P2, mesh, psi1, psi2  = refine_mesh_pureDef(mesh,W,S,old_q,eps,h0,h1,interface_factor,H1)
old_q                           = interpolate(initial, W)

output_mesh(mesh,filepath)
output_solution(old_q,0,mesh,P2,S,psi1,psi2,filepath,0)
