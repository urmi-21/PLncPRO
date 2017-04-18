#!/bin/sh

# BUILD THE MAN PAGES

INPUTDIR=input
OUTPUTDIR=man

HEADERFILE=${INPUTDIR}/ESTate.header.man
FOOTERFILE=${INPUTDIR}/ESTate.footer.man

mkdir $OUTPUTDIR
mkdir ${OUTPUTDIR}/man1
mkdir ${OUTPUTDIR}/man7

ls -1 ${INPUTDIR}/*.[17] | while read FILENAME
do
    MANPAGENAME=`basename $FILENAME`
    echo Building man page documentation for $MANPAGENAME
    MANPAGESECT=`echo $MANPAGENAME | cut -d. -f2`
    cat $HEADERFILE $FILENAME $FOOTERFILE > \
    ${OUTPUTDIR}/man${MANPAGESECT}/$MANPAGENAME
done

echo ---------------------------------------------------------
echo The man pages have been generated.
echo To view them under bash:
echo 
echo export MANPATH=\"`pwd`/man/:\$MANPATH\"
echo 
echo To view them under tcsh:
echo
echo setenv MANPATH \"`pwd`/man/:\$MANPATH\"
echo
echo then type \'man ESTate\' to start view the first manpage.
echo ---------------------------------------------------------
