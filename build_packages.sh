#! /bin/sh
python3 setup.py --command-packages=stdeb.command bdist_deb
cp deb_dist/*.deb packages/
python3 setup.py bdist_rpm
cp dist/*.rpm packages/
