include ./Makedefs

CFLAGS = $(CFDEF)
CLIBS = -lm

all:	pnp_funnel

pnp_funnel:	pnp_funnel_1.o pnp_funnel_2.o pnp_funnel_main.o
	$(CC) -o $@ pnp_funnel_1.o pnp_funnel_2.o pnp_funnel_main.o $(CLIBS)
	strip $@

clean:
	rm -f *.o core pnp_funnel

generic:
	rm -f Makedefs
	ln -s Makedefs.generic Makedefs

alpha:
	rm -f Makedefs
	ln -s Makedefs.alpha Makedefs

hp:
	rm -f Makedefs
	ln -s Makedefs.hp Makedefs

ibm:
	rm -f Makedefs
	ln -s Makedefs.ibm Makedefs

iris:
	rm -f Makedefs
	ln -s Makedefs.iris Makedefs

linux:
	rm -f Makedefs
	ln -s Makedefs.linux Makedefs

sun:
	rm -f Makedefs
	ln -s Makedefs.sun Makedefs

solaris:
	rm -f Makedefs
	ln -s Makedefs.solaris Makedefs
