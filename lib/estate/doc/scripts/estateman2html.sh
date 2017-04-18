#!/bin/sh

# CONVERSION OF MAN PAGES TO HTML
# CHANGES OUTPUT FROM man2html A BIT.

# YOU NEED TO HAVE man2html IN YOUR PATH FOR THIS TO WORK.

STARTDIR=`pwd`
RESULTDIR="html"

INEXPR1='http://localhost/cgi-bin/man/man2html?[1-9]+'
OUTEXPR1='http://www.hgmp.mrc.ac.uk/~gslater/estateman/'

which man2html > /dev/null
if [ $? -eq 0 ]
then
    echo man2html program located ... proceeding
else
    echo No man2html program found - will not make html documentation
    exit 1
fi

mkdir $RESULTDIR

ls -1 man/man[17]/[a-zA-Z]*.[17] | while read MANPAGE
do
     cd `dirname $MANPAGE`
     MANPAGENAME=`basename $MANPAGE`
     MANPAGEHTML=`echo $MANPAGENAME | cut -d'.' -f1`.html
     echo Building html documentation for $MANPAGEHTML
     man2html $MANPAGENAME | tail +2 | sed \
     -e 's/<BODY>/<BODY BGCOLOR=\"ffffff\"><TT>/g' \
     -e 's@</BODY>@</TT></BODY>@g' \
     -e "s@${INEXPR1}\(.*\)\"@${OUTEXPR1}\1.html\"@g" \
     > ${STARTDIR}/${RESULTDIR}/$MANPAGEHTML
     cd $STARTDIR
done

