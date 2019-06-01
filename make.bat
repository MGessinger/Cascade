gcc -O2 -Wall -pedantic -Wextra -o monodromy monodromy.c acb_ode.c -lflint -larb -lgmp
monodromy.exe data/odetest.txt 10
pause
