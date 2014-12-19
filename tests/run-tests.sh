#!/bin/sh
mtree="$1"
testdir="$(dirname $0)"
"$mtree" "$testdir"/binary.nex || true