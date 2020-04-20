#! /bin/sh

# build man
argparse-manpage --pyfile ./bin/metarun  --function get_parser --project-name metadynamic --author R.Plasson --author-email raphael.plasson@univ-avignon.fr  >> docs/man/metarun.1

# build doc
pydoctor --html-output=./docs/apidocs/ metadynamic

# build packages
python3 setup.py --command-packages=stdeb.command bdist_deb
cp deb_dist/*.deb packages/
python3 setup.py bdist_rpm
cp dist/*.rpm packages/
mv *.tar.gz packages/
