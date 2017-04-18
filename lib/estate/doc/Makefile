# ------------------------------------------------------------------
# Makefile for ESTate documentation.
# ------------------------------------------------------------------

# ------------------------------------------------------------------

all:	 man

clean:	cleanman cleanhtml

cleanman:
	rm -rf man

cleanhtml:
	rm -rf html

man:	cleanman
	scripts/buildmanpages.sh;

html:	cleanhtml
	scripts/estateman2html.sh

# ------------------------------------------------------------------

