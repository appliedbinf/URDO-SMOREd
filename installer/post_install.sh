#!/bin/bash

${1}/bin/pip install --upgrade --trusted-host pip.appliedbinf.com --extra-index-url https://cdc:smoredInstaller@pip.appliedbinf.com SMOREd
echo "alias smored=${1}/bin/smored" >> ~/.bashrc

echo ""
echo "================================================="
echo "Installation complete, please reload your shell or run"
echo "the below command to add the smored command"
echo "alias smored=$1/bin/smored"
