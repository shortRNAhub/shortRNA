## @Author: Deepak Tanwar
## @Date:   20180412

##---------##
# This script can find README in markdown format and add TOC to them!
##---------##

#!/usr/bin/env bash

for i in `find . -name "README.md"`
do
    doctoc --github --notitle $i
    if grep -qFi "zenhub" $i; then
	next;
    else
	sed -i '1i<a href="https://zenhub.com"><img src="https://raw.githubusercontent.com/ZenHubIO/support/master/zenhub-badge.png"></a>' $i
    fi
done
