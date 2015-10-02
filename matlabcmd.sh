#!/bin/sh
# This version shows additional error output
script=`basename "$1" .m` # this will return $1 unless it ends with .m
if [ "$script" = "$1" ]; then
    cmdterm=,
else
    cmdterm=\;
fi;
shift # remove the first argument

args="$@"
argsq=`echo $args | tr '"' "'"`
matlab -nodisplay -r "disp('BEGIN>>'); try, $script $argsq $cmdterm catch me, fprintf(2,getReport(me)); exit(1); end, exit(0)" -nosplash | sed -e '1,/^BEGIN>>$/ d'
