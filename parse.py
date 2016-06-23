#!/usr/bin/env python

import ijson

fn = 'test-registry-human.json'

with open(fn) as f:
    for item in ijson.items(f, "item"):
        print item
        break
