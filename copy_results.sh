#!/bin/sh

cd /home/tanburn/Quantum-Factorization-By-Minimization/

cp -r results/20x20/ ../public_html/
cp -r results/30x30/ ../public_html/
cp -r results/40x40/ ../public_html/
cp -r results/50x50/ ../public_html/
cp -r results/60x60/ ../public_html/
cp -r results/70x70/ ../public_html/
cp -r results/80x80/ ../public_html/
cp -r results/90x90/ ../public_html/
cp -r results/100x100/ ../public_html/
cp -r results/110x110/ ../public_html/
cp -r results/120x120/ ../public_html/
cp -r results/130x130/ ../public_html/
cp -r results/140x140/ ../public_html/
cp -r results/150x150/ ../public_html/
cp -r results/160x160/ ../public_html/
cp -r results/170x170/ ../public_html/
cp -r results/180x180/ ../public_html/
cp -r results/190x190/ ../public_html/
cp -r results/200x200/ ../public_html/
cp -r results/210x210/ ../public_html/
cp -r results/220x220/ ../public_html/
cp -r results/230x230/ ../public_html/
cp -r results/240x240/ ../public_html/
cp -r results/250x250/ ../public_html/
./print_results.sh > ../public_html/results_summary.txt
