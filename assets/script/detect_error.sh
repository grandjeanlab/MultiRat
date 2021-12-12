#!/bin/bash

log_path="/scratch/m/mchakrav/gjoanes/multirat/test/export/log/confound/"
find $log_path/*.log | while read line
do
RESULT=`grep 'Error' $line`
if [ -n "$RESULT" ]; then
echo $line
fi
done

exit

