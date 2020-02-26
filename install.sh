#!/bin/bash

python3 setup.py --cython install --record installation.txt >&1 | tee record.txt
