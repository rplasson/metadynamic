#! /bin/bash

# test
mpirun -c 4 --oversubscribe  pytest-3 --with-mpi  docs/tests
if [[ $? -ne 0 ]] ; then
    exit 1
fi


# build man
argparse-manpage --pyfile ./bin/metarun  --function get_parser --project-name metadynamic --author R.Plasson --author-email raphael.plasson@univ-avignon.fr  >> docs/man/metarun.1
sed -i 's/argparse-manpage/metarun/g' docs/man/metarun.1  # Couldn't figure out correct argparse-manpage parameters...

# build doc
pydoctor --html-output=./docs/apidocs/ metadynamic

# build packages
rm -r deb_dist/*
python3 setup.py --command-packages=stdeb.command bdist_deb
cp deb_dist/*.deb packages/
python3 setup.py bdist_rpm
cp dist/*.rpm packages/
mv *.tar.gz packages/
