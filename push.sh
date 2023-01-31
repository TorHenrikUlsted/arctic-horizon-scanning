#!/bin/sh
cd "C:\Users\tor_h\OneDrive - NTNU\Master thesis\MasterCode\MasterThesis-NaturalizingPlants-Arctic"
git checkout dev
git add .
git commit -am "made changes"
git push GitHub dev

cd "C:\Users\tor_h\OneDrive - NTNU\Master thesis\MasterCode\MasterThesis-NaturalizingPlants-Arctic"
<<<<<<< HEAD
git checkout dev
=======
git checkout master
git merge dev
>>>>>>> dev
git add .
git commit -am "made changes"
git push earthengine master

echo Press Enter...
read