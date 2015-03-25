#!/bin/bash

grep "rsem.genes.results" $1 | cut -f 4 | tr "::" "\t" | cut -f 1
