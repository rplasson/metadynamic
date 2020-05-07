#! /bin/bash

# test
mpirun -c 4 --oversubscribe  pytest-3 --with-mpi  docs/tests
if [[ $? -ne 0 ]] ; then
    exit 1
fi


# build doc
pydoctor --html-output=./docs/apidocs/ --project-name=metadynamic metadynamic
#cd docs/sphinx
#make clean
#make html
#make epub
#make latexpdf
#cd -

# build packages
rm -r deb_dist/*
python3 setup.py --command-packages=stdeb.command bdist_deb
cp deb_dist/*.deb packages/
python3 setup.py bdist_rpm
cp dist/*.rpm packages/
mv *.tar.gz packages/
