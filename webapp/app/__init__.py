#!/usr/bin/env python
# Copyright (C) 2015,2016 Mohammad Alanjary
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
# Funding by the German Centre for Infection Research (DZIF)
#
# This file is part of ARTS
# ARTS is free software. you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version
#
# License: You should have received a copy of the GNU General Public License v3 with ARTS
# A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import os
from flask import Flask
from flask_mail import Mail

app = Flask(__name__)

parentdir = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-2])
if os.path.exists(os.path.join(parentdir,"config","active_config.py")):
    app.config.from_pyfile(os.path.join(parentdir,"config","active_config.py"))
else:
    app.config.from_pyfile(os.path.join(parentdir,"config","webapp.py"))
mail = Mail(app)

from app import views