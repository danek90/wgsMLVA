#!/bin/sh

DATE=`date +%Y-%m-%d`

git add -A;
git commit -m "update: $DATE";
git push https://github.com/danek90/fragments