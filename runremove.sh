#!/bin/bash

tr '\n' '\0' < installation.txt | xargs -0 rm -f --
rm -r dist
rm -r build
rm -r owcsimpy.egg-info/
