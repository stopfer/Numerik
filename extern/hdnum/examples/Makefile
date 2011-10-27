CC          = g++
HDNUMPATH   = ..
#CCOPTS      = -I$(HDNUMPATH) -g -O0
CCOPTS      = -I$(HDNUMPATH) -O3
LIOPTSGMP   = -lm -L/opt/local/lib -lgmpxx -lgmp
LIOPTS      = -lm

hallowelt_MAIN = hallowelt
hallowelt_SOURCES =

zahlen_MAIN = zahlen
zahlen_SOURCES =

string_MAIN = string
string_SOURCES =

mehrzahlen_MAIN = mehrzahlen
mehrzahlen_SOURCES =

eingabe_MAIN = eingabe
eingabe_SOURCES =

pendelwhile_MAIN = pendelwhile
pendelwhile_SOURCES =

pendel_MAIN = pendel
pendel_SOURCES =

pendel_MAIN = pendelmittimer
pendel_SOURCES =

pendelnumerisch_MAIN = pendelnumerisch
pendelnumerisch_SOURCES =

pendelmitfunktion_MAIN = pendelmitfunktion
pendelmitfunktion_SOURCES =

pendelmitfunktionstemplate_MAIN = pendelmitfunktionstemplate
pendelmitfunktionstemplate_SOURCES =

whitespace_MAIN = whitespace
whitespace_SOURCES =

funktion_MAIN = funktion
funktion_SOURCES =

precision_MAIN = precision
precision_SOURCES =

hallohdnum_MAIN = hallohdnum
hallohdnum_SOURCES =

fibonacci_MAIN = fibonacci
fibonacci_SOURCES =

vektoren_MAIN = vektoren
vektoren_SOURCES =

matrizen_MAIN = matrizen
matrizen_SOURCES =

lr_MAIN = lr
lr_SOURCES =

modelproblem_MAIN = modelproblem
modelproblem_SOURCES =

ordertest_MAIN = ordertest
ordertest_SOURCES =

hodgkinhuxley_MAIN = hodgkinhuxley
hodgkinhuxley_SOURCES =

nbody_MAIN = nbody
nbody_SOURCES =

vanderpol_MAIN = vanderpol
vanderpol_SOURCES =

wurzel_MAIN = wurzel
wurzel_SOURCES =

iemodelproblem_MAIN = iemodelproblem
iemodelproblem_SOURCES =

ecke_MAIN = ecke
ecke_SOURCES =

all : nogmp gmp

nogmp : $(hallowelt_MAIN) $(zahlen_MAIN) $(string_MAIN)  $(pendelwhile_MAIN) $(pendel_MAIN) $(pendelmittimer_MAIN) $(pendelnumerisch_MAIN) $(pendelmitfunktion_MAIN) $(pendelmitfunktionstemplate_MAIN) $(whitespace_MAIN) $(funktion_MAIN) $(precision_MAIN) $(hallohdnum_MAIN) $(vektoren_MAIN) $(matrizen_MAIN) $(modelproblem_MAIN) $(ordertest_MAIN) $(hodgkinhuxley_MAIN) $(nbody_MAIN)  $(vanderpol_MAIN) $(wurzel_MAIN) $(iemodelproblem_MAIN) $(fibonacci_MAIN)

gmp : $(mehrzahlen_MAIN) $(eingabe_MAIN) $(lr_MAIN)

$(hallowelt_MAIN): $(hallowelt_MAIN).o $(hallowelt_SOURCES) Makefile
	$(CC) -o $(hallowelt_MAIN) $(hallowelt_MAIN).o $(hallowelt_SOURCES) $(LIOPTS)

$(zahlen_MAIN): $(zahlen_MAIN).o $(zahlen_SOURCES) Makefile
	$(CC) -o $(zahlen_MAIN) $(zahlen_MAIN).o $(zahlen_SOURCES) $(LIOPTS)

$(string_MAIN): $(string_MAIN).o $(string_SOURCES) Makefile
	$(CC) -o $(string_MAIN) $(string_MAIN).o $(string_SOURCES) $(LIOPTS)

$(pendelwhile_MAIN): $(pendelwhile_MAIN).o $(pendelwhile_SOURCES) Makefile
	$(CC) -o $(pendelwhile_MAIN) $(pendelwhile_MAIN).o $(pendelwhile_SOURCES) $(LIOPTS)

$(pendel_MAIN): $(pendel_MAIN).o $(pendel_SOURCES) Makefile
	$(CC) -o $(pendel_MAIN) $(pendel_MAIN).o $(pendel_SOURCES) $(LIOPTS)

$(pendelmittimer_MAIN): $(pendelmittimer_MAIN).o $(pendelmittimer_SOURCES) Makefile
	$(CC) -o $(pendelmittimer_MAIN) $(pendelmittimer_MAIN).o $(pendelmittimer_SOURCES) $(LIOPTS)

$(pendelnumerisch_MAIN): $(pendelnumerisch_MAIN).o $(pendelnumerisch_SOURCES) Makefile
	$(CC) -o $(pendelnumerisch_MAIN) $(pendelnumerisch_MAIN).o $(pendelnumerisch_SOURCES) $(LIOPTS)

$(pendelmitfunktion_MAIN): $(pendelmitfunktion_MAIN).o $(pendelmitfunktion_SOURCES) Makefile
	$(CC) -o $(pendelmitfunktion_MAIN) $(pendelmitfunktion_MAIN).o $(pendelmitfunktion_SOURCES) $(LIOPTS)

$(pendelmitfunktionstemplate_MAIN): $(pendelmitfunktionstemplate_MAIN).o $(pendelmitfunktionstemplate_SOURCES) Makefile
	$(CC) -o $(pendelmitfunktionstemplate_MAIN) $(pendelmitfunktionstemplate_MAIN).o $(pendelmitfunktionstemplate_SOURCES) $(LIOPTS)

$(whitespace_MAIN): $(whitespace_MAIN).o $(whitespace_SOURCES) Makefile
	$(CC) -o $(whitespace_MAIN) $(whitespace_MAIN).o $(whitespace_SOURCES) $(LIOPTS)

$(funktion_MAIN): $(funktion_MAIN).o $(funktion_SOURCES) Makefile
	$(CC) -o $(funktion_MAIN) $(funktion_MAIN).o $(funktion_SOURCES) $(LIOPTS)

$(precision_MAIN): $(precision_MAIN).o $(precision_SOURCES) Makefile
	$(CC) -o $(precision_MAIN) $(precision_MAIN).o $(precision_SOURCES) $(LIOPTS)

$(hallohdnum_MAIN): $(hallohdnum_MAIN).o $(hallohdnum_SOURCES) Makefile
	$(CC) -o $(hallohdnum_MAIN) $(hallohdnum_MAIN).o $(hallohdnum_SOURCES) $(LIOPTS)

$(vektoren_MAIN): $(vektoren_MAIN).o $(vektoren_SOURCES) Makefile
	$(CC) -o $(vektoren_MAIN) $(vektoren_MAIN).o $(vektoren_SOURCES) $(LIOPTS)

$(matrizen_MAIN): $(matrizen_MAIN).o $(matrizen_SOURCES) Makefile
	$(CC) -o $(matrizen_MAIN) $(matrizen_MAIN).o $(matrizen_SOURCES) $(LIOPTS)

$(modelproblem_MAIN): $(modelproblem_MAIN).o $(modelproblem_SOURCES) Makefile
	$(CC) -o $(modelproblem_MAIN) $(modelproblem_MAIN).o $(modelproblem_SOURCES) $(LIOPTS)

$(fibonacci_MAIN): $(fibonacci_MAIN).o $(fibonacci_SOURCES) Makefile
	$(CC) -o $(fibonacci_MAIN) $(fibonacci_MAIN).o $(fibonacci_SOURCES) $(LIOPTS)

$(ordertest_MAIN): $(ordertest_MAIN).o $(ordertest_SOURCES) Makefile
	$(CC) -o $(ordertest_MAIN) $(ordertest_MAIN).o $(ordertest_SOURCES) $(LIOPTS)

$(hodgkinhuxley_MAIN): $(hodgkinhuxley_MAIN).o $(hodgkinhuxley_SOURCES) Makefile
	$(CC) -o $(hodgkinhuxley_MAIN) $(hodgkinhuxley_MAIN).o $(hodgkinhuxley_SOURCES) $(LIOPTS)

$(nbody_MAIN): $(nbody_MAIN).o $(nbody_SOURCES) Makefile
	$(CC) -o $(nbody_MAIN) $(nbody_MAIN).o $(nbody_SOURCES) $(LIOPTS)

$(vanderpol_MAIN): $(vanderpol_MAIN).o $(vanderpol_SOURCES) Makefile
	$(CC) -o $(vanderpol_MAIN) $(vanderpol_MAIN).o $(vanderpol_SOURCES) $(LIOPTS)

$(wurzel_MAIN): $(wurzel_MAIN).o $(wurzel_SOURCES) Makefile
	$(CC) -o $(wurzel_MAIN) $(wurzel_MAIN).o $(wurzel_SOURCES) $(LIOPTS)

$(iemodelproblem_MAIN): $(iemodelproblem_MAIN).o $(iemodelproblem_SOURCES) Makefile
	$(CC) -o $(iemodelproblem_MAIN) $(iemodelproblem_MAIN).o $(iemodelproblem_SOURCES) $(LIOPTS)

$(ecke_MAIN): $(ecke_MAIN).o $(ecke_SOURCES) Makefile
	$(CC) -o $(ecke_MAIN) $(ecke_MAIN).o $(ecke_SOURCES) $(LIOPTS)


# Applications with GNU multiprecision library

$(eingabe_MAIN): $(eingabe_MAIN).o $(eingabe_SOURCES) Makefile
	$(CC) -o $(eingabe_MAIN) $(eingabe_MAIN).o $(eingabe_SOURCES) $(LIOPTSGMP)

$(mehrzahlen_MAIN): $(mehrzahlen_MAIN).o $(mehrzahlen_SOURCES) Makefile
	$(CC) -o $(mehrzahlen_MAIN) $(mehrzahlen_MAIN).o $(mehrzahlen_SOURCES) $(LIOPTSGMP)

$(lr_MAIN): $(lr_MAIN).o $(lr_SOURCES) Makefile
	$(CC) -o $(lr_MAIN) $(lr_MAIN).o $(lr_SOURCES) $(LIOPTSGMP)

.cc.o:
	$(CC) $(CCOPTS) -c $<

clean:
	rm -f *.o *.dat $(hallowelt_MAIN) $(zahlen_MAIN) $(string_MAIN) $(mehrzahlen_MAIN) $(eingabe_MAIN) $(pendelwhile_MAIN) $(pendel_MAIN) $(pendelnumerisch_MAIN) $(pendelmitfunktion_MAIN) $(pendelmitfunktionstemplate_MAIN) $(whitespace_MAIN) $(funktion_MAIN) $(precision_MAIN) $(hallohdnum_MAIN) $(vektoren_MAIN) $(matrizen_MAIN) $(modelproblem_MAIN) $(ordertest_MAIN) $(hodgkinhuxley_MAIN) $(nbody_MAIN) $(lr_MAIN) $(vanderpol_MAIN) $(wurzel_MAIN) $(iemodelproblem_MAIN)

