#!/bin/sh
exe="$1"
t="$(dirname $0)"
passed=0
failed=0
function runcheck
{
    if python "${t}/check-lnL.py" "${exe}" "$1" "$2" "$3"
    then
        passed=$(expr $passed + 1)
    else
        failed=$(expr $failed + 1)
    fi
}

runcheck "${t}/jc.nex" "${t}/jc.ini" "-76.66761"
if test $failed -gt 0
then
    echo "Failed ${failed} / " $(expr $passed + $failed) " tests."
else
    echo "Passed all $passed tests."
fi
