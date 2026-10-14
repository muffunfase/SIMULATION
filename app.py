import sympy as sp
from sympy import sin, cos
import numpy as np
import math

import csv
import pandas as pandas
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

import os

import gmsh
if gmsh.isInitialized():
    gmsh.finalize()          # kill any leftover model from a previous run
gmsh.initialize()
gmsh.clear()                 # empty model, just to be sure
if os.path.exists("TEST.msh"):
    os.remove("TEST.msh")  

import sqlite3

from flask import g
from flask import Flask, render_template, request, jsonify, redirect, url_for, request

if os.path.exists("app.db"):
    os.remove("app.db")

app = Flask(__name__)

dbfile = 'app.db'

def get_db():
    connection = g.get('db', 'null')
    if connection == 'null':
        g.db = sqlite3.connect(dbfile)
        g.db.row_factory=sqlite3.Row
        return g.db
    else:
        return connection

def create_table():
    connection = get_db()
    sql=connection.cursor()
    connection.commit()
    sql.execute('''CREATE TABLE IF NOT EXISTS ic (
    "id" integer primary key autoincrement,
    'name' TEXT NOT NULL,
    "density" REAL DEFAULT 998.2,
    "dynVic" Real DEFAULT 1.002e-3,
    'radius' REAL DEFAULT 0.005,
    'length' REAL DEFAULT 0.1,
    'inVel' REAL DEFAULT 0.005,
    'inProf' TEXT DEFAULT 'uniform',
    'outP' REAL DEFAULT 0,
    'Wall_CON' TEXT DEFAULT 'no_slip',

    'MeshSz' REAL DEFAULT 0.001,
    'VelOrd' INTEGER DEFAULT 2,
    'POrd' INTEGER DEFAULT 1,

    'model' TEXT DEFAULT 'navier stokes',
    'initGuess' TEXT DEFAULT 'stokes',
    'MaxIT' INTEGER DEFAULT 50,
    'Tol' REAL DEFAULT 1e-8,

    'timedep' INTEGER DEFAULT 0,
    'dt' REAL DEFAULT NULL,
    't_f' REAL DEFAULT NULL,
    'vel_i' REAL DEFAULT 0.0

    )''')

def table_insert(name, density):
    with app.app_context():
        connection = get_db()
        sql = connection.cursor()
        insertion = [name, density]
        sql.execute('''
        INSERT INTO ic (name, density) VALUES (?,?)
        ''', list(insertion))
        connection.commit()

with app.app_context():
    create_table()


def get_db_val():
    with app.app_context():
        connection = get_db()
        sql=connection.cursor()
        connection.commit()
        datum = sql.execute('SELECT * FROM ic')
        for item in datum:
            for col in item:
                print(col)

with app.app_context():
    get_db_val()

@app.route('/')
def home():
    return render_template('index.html')

table_insert('TEST', 123)
get_db_val()



## MESHING

### Initializations
R   = 0.005       # pipe radius, m    (must match the mesh)
L   = 0.05        # pipe length, m    (must match the mesh)
h = 0.0015
rho = 998.2       # water density, kg/m³
mu  = 1.002e-3    # water dynamic viscosity, Pa·s
U   = 0.005       # mean inlet velocity, m/s
MaxIT = 50
Tol   = 1e-8
#


gmsh.initialize()
gmsh.model.add('TEST')
occ = gmsh.model.occ
inside = occ.addCylinder(0,0,0, 0,0,L, R )
occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", h)   # element size, in the same units as the geometry
gmsh.model.mesh.generate(3)
gmsh.write('TEST.msh')
gmsh.fltk.run()


######################################################################################################
from skfem import *
from skfem.helpers import grad, div, dot, ddot, mul
from scipy.sparse import bmat
m = MeshTet.load("TEST.msh")


m = m.with_boundaries({
    "inlet":  lambda x: np.isclose(x[2], 0),
    "outlet": lambda x: np.isclose(x[2], L),
    "wall":   lambda x: ~np.isclose(x[2], 0) & ~np.isclose(x[2], L),
})



ub = Basis(m, ElementVectorH1(ElementTetP2()), intorder=4) ##Velocity
pb = Basis(m, ElementTetP1(), intorder=4) ## Pressure
nu_, np_ = ub.N, pb.N

@BilinearForm
def visc(u, v, w):
    return mu * ddot(grad(u), grad(v))

@BilinearForm
def divf(u, q, w):
    return -div(u) * q

A = visc.assemble(ub)
B = divf.assemble(ub, pb)
K = bmat([[A, B.T], [B, None]], "csr")

x = np.zeros(nu_ + np_)
D_in = ub.get_dofs('inlet')
D_wall = ub.get_dofs('wall')

iz = D_in.all("u^3")
xd, yd = ub.doflocs[0, iz], ub.doflocs[1, iz]
x[iz] = 2*U*(1 - (xd**2 + yd**2)/R**2)

x[D_wall.all()] = 0.0
Dfix = np.concatenate([D_in.all(), D_wall.all()])

x = solve(*condense(K, x=x, D=Dfix))
u, p = x[:nu_], x[nu_:]

uz_basis = ub.split_bases()[2]
uz = u[ub.split_indices()[2]]

r = np.array([0, 0.25 , 0.5, 0.75, 0.95]) * R
pts = np.vstack([r, 0*r, 0*r + 0.9*L])
print(uz_basis.probes(pts) @ uz)

pz = pb.probes(np.array([[0, 0, 0.6*L], [0, 0, 0.9*L]]).T) @ p
print((pz[1] - pz[0]) / (0.3*L), "vs exact", -8*mu*U/R**2) 
print('Error:', abs((pz[1] - pz[0]) / (0.3*L) - (-8*mu*U/R**2)))



import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(8, 4))
gs = gridspec.GridSpec(3,4)

axlSIMULATION = fig.add_subplot(gs[0:2, 0:2])
axlSIMULATION.set_title('SIMULATION')
axlHEATMAP = fig.add_subplot(gs[0, 2])
axlHEATMAP.set_title('HEAT MAP')
axlRADIATION = fig.add_subplot(gs[0, 3])
axlRADIATION.set_title('RADIATION HEAT MODEL')
axlBLUEPRINT = fig.add_subplot(gs[1:2, 2:4])
axlBLUEPRINT.set_title('BLUEPRINT')

plt.show()

### Goals:
# SIM MUST HAVE:
# Simulation of Navier-Stokes Maxwell equations CFD in the MPD thruster
# Heat Map
# Radiation Heat Model
# Realistic Simulation
# Blueprint
# Legend showing Efficiency, velocity, acceleration, and thrust