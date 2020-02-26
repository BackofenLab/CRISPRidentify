#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 indexname"
  echo "migrate a virtual index from Sparc/Alpha- to Intel-platforms and vice versa."
  exit 1
fi

EXTENSIONS="bck llv suf sds ssp skp"
TMPFILE=/tmp/.tmp
indexname=$1

if test -f ${indexname}.prj
then
  echo "# migrate the indexfile ${indexname}.{$EXTENSIONS}"
else
  echo "$0: index \"${indexname}\" does not exist"
  exit 1
fi

for ext in $EXTENSIONS
do 
  if test -f ${indexname}.${ext}
  then
    echo "Converting ${indexname}.${ext}"
    vendian 4 ${indexname}.${ext} > $TMPFILE
    if test $? -ne 0
    then
      echo "failure: vendian 4 ${indexname}.${ext}"
      exit 1
    fi
    mv $TMPFILE ${indexname}.${ext}
  fi
done

echo "Done."
