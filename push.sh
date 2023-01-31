#!/bin/sh
cd "C:\Users\tor_h\OneDrive - NTNU\Master thesis\MasterCode\MasterThesis-NaturalizingPlants-Arctic"
git checkout dev
git add .
git commit -am "made changes"
git push GitHub dev

git checkout dmaster
git add .
git commit -am "made changes"
git push earthengine master
echo Press Enter...
read