   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%                       ------------                        %%
   %%                       ESTate 0.5.0                        %%
   %%                       ------------                        %%
   %%     >>> Expressed Sequence Tag Analysis Tools Etc <<<     %%
   %%                (c) Guy St.C. Slater 1996-1999.            %%
   %%                                                           %%
   %%              Human Genome Mapping Project RC,             %%
   %%                   Hinxton, Cambridge, UK.                 %%
   %%                                                           %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q. What is ESTate ?

A. ESTate stands for Expressed Sequence Tag Analysis Tools Etc.

   It contains a bunch of programs for EST clustering,
   database searching, and sequence translation, and so on.

-//-

Q. Where is the documentation ?

A. See the man pages.  man ESTate is a good place to start.

   To view the man pages without installing them on your system,

   For csh,
   % setenv MANPATH "${ESTATE}/doc/man:$MANPATH"

   For bash,
   % export MANPATH="${ESTATE}/doc/man:$MANPATH"

   where $ESTATE is the direcory where you unpacked the distribution.

-//-

Q. How do I install ESTate ?

   Just type 'make' in the top level directory.
   This should build the programs and documenation.

   You may wish to edit the Makefile in the src directory
   to customise compiler optimisations etc.

   ESTate has been developed on IRIX and Linux systems,
   but has also been compiled on Solaris and AIX.

   Please let me know if you have problems compiling.

-//-

Q.  Under what terms is the source code available ?

A.  The GNU public license.  See the file COPYING for details,
    or http://www.fsf.org/copyleft/gpl.html for more information.

-//-

Q.  Problems ?  Feedback ?

A.  All feedback is appreciated.
    Please mail me at gslater@hgmp.mrc.ac.uk
--
