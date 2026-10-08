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
