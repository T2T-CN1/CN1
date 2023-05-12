#!/usr/bin/bash

git add $1
read -p "Enter your comments for this change: " comm
git commit -a -m "$comm"
# some are master, some are main
# git push origin main
# git push origin master
git push
