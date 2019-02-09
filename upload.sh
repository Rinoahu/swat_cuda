#!/bin/bash
#cd MCL_lite
#rm -rf __pycache__
#rm module/tests/tmp*.npz* mcl_sparse.pyc
git config --global user.email xiaohu@iastate.edu
git config --global user.name Rinoahu


git remote rm origin

git add -A .
git commit -m 'add cpu version of parallelized swat to verify my parallization of swat algorithm'
git remote add origin https://github.com/Rinoahu/swat_cuda

git pull origin master
git push origin master

git checkout master
