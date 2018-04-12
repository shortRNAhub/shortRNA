## @Author: Deepak Tanwar
## @Date:   20180412

##---------##
# This script can find README in markdown format and add TOC to them!
##---------##

#!/usr/bin/env bash

for i in `find . -name "README.md"`
do
    doctoc --github --notitle $i
done
