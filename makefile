
#    Authors :
#    KAOUCH Abdelssamad
#    AIT KHELIFA Tanina


all: ver_double ver_mpfr

ver_double: ver_double.c
	gcc -o ver_double ver_double.c -Wall -lmpfr -lgmp -lm

ver_mpfr: ver_mpfr.c
	gcc -o ver_mpfr ver_mpfr.c -Wall -lmpfr -lgmp -lm

clean:
	-rm -f ver_mpfr.o
	-rm -f ver_mpfr
	-rm -f ver_double.o
	-rm -f ver_double